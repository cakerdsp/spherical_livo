# Spherical LIVO — ROS1

基于官方 FAST-LIVO2 的多相机球面直接法实现。基准提交：
`0d2c0346107b75b59934975adec9a6eeeb913c64`。

**当前交付状态：源码与基本静态检查；未编译、未运行测试、未回放数据。**
本机没有 ROS1 编译条件。本文不将实现状态表述为精度或实时性能已验证。

## 范围

- 保留 FAST-LIVO2 的雷达预处理、IMU 传播／去畸变和体素平面地图。
- 多相机共用一套球面直接法、一张视觉地图，历史参考可以来自其他相机。
- 原始像素仅用于采集数据和 LUT 索引。采样节点、插值核、梯度、warp 和视觉残差均使用三维单位方向。
- 雷达局部平面给每个采样方向提供三维支持。跨深度断层、遮挡证据、无可靠平面支持和饱和支持均拒绝使用。
- 无光流前端、虚拟针孔图像、在线外参／时间标定、方向性更新或高级视觉地图管理。
- 复杂光度模型与专门弱观测处理见 [研究 TODO](docs/research_todo.md)。

## 代码入口

| 文件 | 责任 |
| --- | --- |
| `src/spherical.cpp` | 相机射线 LUT、方向查询、球面滤波插值与一致梯度 |
| `src/spherical_vio.cpp` | 平面支持、整块 warp、带宽匹配、曝光联合估计、视觉更新与共享参考 |
| `src/main.cpp` | ROS1 输入、按时间处理相机事件、雷达／IMU 接入、结果输出 |
| `src/IMU_Processing.cpp`、`src/voxel_map.cpp`、`src/preprocess.cpp` | 来自 FAST-LIVO2 的传感器主干 |

详细模型、近似与固定常数记录在 [方法说明](docs/method.md)。

## ROS1 环境

完整命令见 **[ROS1 编译与部署教程](docs/deployment_ros1_zh.md)**，包含 Ubuntu 20.04 / Noetic、SDK2 与驱动依赖、三种数据集回放和参数依据。必须复制这份已修改源码；重新克隆上游 FAST-LIVO2 不包含本实现。

在完成教程中的依赖安装后：

```bash
source /opt/ros/noetic/setup.bash
source ~/livox_ros1_ws/devel/setup.bash
cd ~/spherical_ws
catkin_make -DCMAKE_BUILD_TYPE=Release -DCATKIN_ENABLE_TESTING=OFF -j2 -l2
source devel/setup.bash
roslaunch spherical_livo mapping_m2dgr.launch num_cameras:=3 output_dir:=$HOME/spherical_results/run01
```

等待 `Spherical LIVO ready`，在另一个已 source 环境的终端播放数据：

```bash
source ~/spherical_ws/devel/setup.bash
rosbag play -r 0.5 /absolute/path/to/sequence.bag
```

所有启动文件默认 `use_sim_time:=false`，普通 `rosbag play` 即可；主循环使用墙上时钟调度，估计器仍使用传感器消息时间戳。需要仿真时间时显式设置 `use_sim_time:=true` 并配合 `rosbag play --clock`。

相同数据重复播放需重新启动节点。输入时间倒退会明确报错，不在已有地图中拼接不同时间段。

## 配置

| 数据集 | 启动文件 | 默认相机 | 参数 / 内参 |
| --- | --- | --- | --- |
| 私有 MID360 + MEI | `mapping_private_mid360.launch` | 左、右两路 | `private_mid360.yaml` / `cameras_private_mid360.yaml` |
| HILTI22 | `mapping_hilti22.launch` | C0、C3、C4 三路，可选五路 | `hilti22.yaml` / `cameras_hilti22.yaml` |
| M2DGR | `mapping_m2dgr.launch` | left、third、fourth 三路，可选五路 | `m2dgr.yaml` / `cameras_m2dgr.yaml` |

所有启动入口默认开启 RViz，共用 `rviz_cfg/spherical_livo.rviz`：固定坐标系 `world`，显示最近 30 秒点云、轨迹和当前 IMU 位姿。无需图形界面时加 `rviz:=false`；关闭 RViz 窗口不会停止估计器。目标机需要安装 `ros-noetic-rviz`。

配置在 `config/` 下。私有与 M2DGR 标定来自 cake_slam；HILTI22 使用重新整理的官方逐相机标定，见 [来源与编号对应](docs/calibration_sources.md)。参数均为未回放验证的默认值。私有 PointCloud2 使用类型 **8**，与旧项目的类型编号不同。 该路径直接移植 cake_slam 的字段解析与时间归一化逻辑，支持 `offset_time/timestamp/time/t` 和 `line/ring`；无有效时间跨度时沿用旧项目按 `scan_rate` 生成近似逐点时间的行为。

自有设备使用：

```bash
roslaunch spherical_livo mapping.launch config:=/absolute/device.yaml cameras:=/absolute/cameras.yaml
```

- `extrin_calib/extrinsic_R,T`：LiDAR → IMU，`p_I = R_IL p_L + t_IL`。
- 每相机的 `Rcl,Pcl`：LiDAR → Camera；`time_offset` 为加到该图像时间戳上的固定秒数。
- `camera_ns` 指向内参节。支持 Pinhole（radtan）、EquidistantCamera（KB4）和 MEI。
- 使用原始标定尺寸、`scale: 1.0`，不做输入去畸变或统一针孔重采样。
- 图像话题可以是原始 Image 或以 `/compressed` 结尾的 CompressedImage。
- 相机按真实时间逐个处理，精确同时间的观测顺序更新同一状态；不把不同时间的图像强行设置为同一时刻。
- IMU/LiDAR 噪声、雷达分辨率和固定标定属于设备配置；图像噪声由当前原图的稳健高频统计估计。

视觉只暴露三个配置：

| 参数 | 默认 | 含义 |
| --- | --- | --- |
| `patch_radius_deg` | 0.7 | 参考球冠角半径 |
| `patch_samples` | 48 | 固定球面采样数量 |
| `max_patches` | 120 | 单次相机事件的 patch 上限 |

不需要设置虚拟焦距、像素 patch 大小、仿射 warp、方向性阈值或按相机逐项配置视觉权重。固定工程常数并未被称为“没有超参数”，其依据和限制见方法说明。

## 输出与基本资源边界

- `/spherical_livo/odometry`：IMU 在 world 下的位姿，包含旋转／平移交叉协方差转换。
- `/spherical_livo/path`、`/spherical_livo/cloud`：轨迹与雷达点云。
- 设置 `output_dir` 后写出 `trajectory.txt`（TUM 格式）和 `visual.csv`。
- `visual.csv` 记录实际视觉更新、跨相机参考数量、深度拒绝、残差和耗时；不能仅凭雷达里程计正常运行判断视觉成功。
- 参考仅存原图局部块，射线 LUT 按相机共享。视觉点最多为 `32 × max_patches`，到容量时释放最早条目；这是基本所有权和内存限制，不含额外质量／生命周期策略。
- LiDAR 可以在相机缺失时继续运行。无待处理图像时保留一帧原生雷达周期的重排缓冲；迟到至已处理时间之前的图像会被丢弃并提示。

## 以后有环境时的验证

仓库包含 `tests/spherical_geometry_test.cpp` 和 `tests/spherical_integration_test.cpp`，覆盖球面导数、平面 warp、饱和支持、负 Z 方向及跨相机共享参考的合成案例。**此次没有执行这些测试。**

```bash
cd ~/spherical_ws
catkin_make run_tests_spherical_livo
catkin_test_results build/test_results
```

实际数据验收仍需记录视觉有效更新与跨相机复用、轨迹误差、逐帧耗时和内存占用。没有在此预设或宣称优于 FAST-LIVO2／cake_slam。

## 来源与许可

基础代码：[hku-mars/FAST-LIVO2](https://github.com/hku-mars/FAST-LIVO2)，GPL-2.0；保留原有源文件许可头及根目录 LICENSE。M2DGR 配置迁移自本地 cake_slam。新视觉前端和 ROS1 集成使用同一许可。
