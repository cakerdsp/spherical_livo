# ROS1 编译与数据集部署

本教程对应本地 `D:\GithubRepo\spherical_livo` 的球面 patch 实现。三套预设已经写入 `config/` 和 `launch/`；参数依据现有标定、原图角分辨率和代码约束选择，尚未经过编译或数据回放验证。

## 1. 部署环境与源码

以下命令面向 **Ubuntu 20.04 + ROS Noetic + GCC 9**，使用 ROS1 catkin。Noetic 的官方目标平台包含 Ubuntu 20.04，其维护期已于 2025 年 5 月结束；这里是为了使用明确的 ROS1 依赖组合，不将它描述为仍受维护的新平台。[ROS REP-3](https://github.com/ros-infrastructure/rep/blob/master/rep-0003.rst#noetic-ninjemys-may-2020---may-2025)

尚未安装 ROS1 的机器先按 [Noetic 安装文档](https://wiki.ros.org/noetic/Installation/Ubuntu) 配置软件源并安装 Noetic。已有 ROS1 环境可直接继续。不要在 source ROS2 Humble 的终端里执行这些命令。

需要转移的是**这份已修改的源码**。目前 `origin` 仍指向官方 FAST-LIVO2，重新克隆它或导出上游 HEAD 都不会得到此次实现。随本次交付提供 `spherical_livo_source.zip`，包含当前源码、配置和文档，不含 `.git`、数据集和构建产物。

在目标 Ubuntu 机器准备好压缩包后：

```bash
mkdir -p ~/spherical_ws/src
unzip ~/Downloads/spherical_livo_source.zip -d ~/spherical_ws/src
```

解压后的文件应为 `~/spherical_ws/src/spherical_livo/package.xml`，不要多套一层同名目录。也可直接复制整个本地 `spherical_livo` 文件夹到 `~/spherical_ws/src/`。

## 2. 安装依赖

下面只是在文档中提供目标机命令；本次交付没有执行安装或编译。

```bash
source /opt/ros/noetic/setup.bash
sudo apt update
sudo apt install -y build-essential cmake git unzip \
  libeigen3-dev libopencv-dev libpcl-dev \
  ros-noetic-roscpp ros-noetic-sensor-msgs ros-noetic-geometry-msgs \
  ros-noetic-nav-msgs ros-noetic-visualization-msgs \
  ros-noetic-pcl-ros ros-noetic-pcl-conversions ros-noetic-tf \
  ros-noetic-cv-bridge ros-noetic-rosbag ros-noetic-roslaunch ros-noetic-rviz
```

没有 GPU/CUDA、Sophus 或 vikit 依赖。预处理保留了 `livox_ros_driver2/CustomMsg` 接口，因此三种数据集都需要安装 **ROS1 版** Livox Driver 2 消息依赖；离线回放时不需要启动驱动节点。

### Livox SDK2 与驱动

若目标机已有编译成功的 ROS1 `livox_ros_driver2` 工作空间，source 它的 `devel/setup.bash`，然后跳到第 3 节即可。否则执行：

```bash
mkdir -p ~/third_party
cd ~/third_party
git clone https://github.com/Livox-SDK/Livox-SDK2.git
cmake -S Livox-SDK2 -B Livox-SDK2/build -DCMAKE_BUILD_TYPE=Release
cmake --build Livox-SDK2/build --parallel 2
sudo cmake --install Livox-SDK2/build
sudo ldconfig

mkdir -p ~/livox_ros1_ws/src
git clone https://github.com/Livox-SDK/livox_ros_driver2.git \
  ~/livox_ros1_ws/src/livox_ros_driver2
cd ~/livox_ros1_ws/src/livox_ros_driver2
cp package_ROS1.xml package.xml
cd ~/livox_ros1_ws
source /opt/ros/noetic/setup.bash
catkin_make -DROS_EDITION=ROS1 -DCMAKE_BUILD_TYPE=Release -j2 -l2
source ~/livox_ros1_ws/devel/setup.bash
```

这里手动执行官方 `build.sh ROS1` 的 ROS1 清单选择和编译步骤，并把驱动放在独立工作空间；官方脚本会删除所在工作空间的 `build/devel/install`，因此不建议在已有算法工作空间中直接重跑它。[SDK2 安装](https://github.com/Livox-SDK/Livox-SDK2#22-instruction-for-ubuntu-2004)、[驱动 build.sh](https://github.com/Livox-SDK/livox_ros_driver2/blob/master/build.sh)

## 3. 编译本项目

```bash
source /opt/ros/noetic/setup.bash
source ~/livox_ros1_ws/devel/setup.bash
cd ~/spherical_ws
catkin_make -DCMAKE_BUILD_TYPE=Release -DCATKIN_ENABLE_TESTING=OFF -j2 -l2
source ~/spherical_ws/devel/setup.bash
```

若使用已有的驱动工作空间，把上面的 `~/livox_ros1_ws` 换成它的路径。首次编译使用两路并行，内存较小可改为 `-j1 -l1`。`CATKIN_ENABLE_TESTING=OFF` 只关闭测试目标；正常估计器仍完整构建。

后续运行的新终端按相同顺序 source ROS1、驱动、算法工作空间；不要再 source 旧 cake_slam 工作空间覆盖包搜索路径。编译成功后，实际启动文件位于包 `spherical_livo` 中。

## 4. 三套直接启动入口

| 数据 | 启动文件 | 默认相机 | 参数 / 内参文件 |
| --- | --- | --- | --- |
| 私有 MID360 + MEI | `mapping_private_mid360.launch` | 左、右，共 2 路 | `private_mid360.yaml` / `cameras_private_mid360.yaml` |
| HILTI22 | `mapping_hilti22.launch` | 物理 C0、C3、C4，共 3 路 | `hilti22.yaml` / `cameras_hilti22.yaml` |
| M2DGR | `mapping_m2dgr.launch` | left、third、fourth，共 3 路 | `m2dgr.yaml` / `cameras_m2dgr.yaml` |

以下三条选择与数据对应的一条。`output_dir` 会自动创建；每次使用不同目录，避免覆盖已有结果。

```bash
roslaunch spherical_livo mapping_private_mid360.launch \
  output_dir:=$HOME/spherical_results/private_run01

roslaunch spherical_livo mapping_hilti22.launch \
  output_dir:=$HOME/spherical_results/hilti22_run01

roslaunch spherical_livo mapping_m2dgr.launch \
  output_dir:=$HOME/spherical_results/m2dgr_run01
```

启动阶段需要为各相机建立原始射线 LUT。等终端出现 `Spherical LIVO ready` 后，在另一个已 source 环境的终端开始回放：

```bash
rosbag play --pause -r 0.5 /absolute/path/to/sequence.bag
```

所有启动入口默认 `use_sim_time:=false`，不要求 `--clock`。主循环按墙上时钟调度，IMU 传播、去畸变和相机同步仍使用消息中的传感器时间戳。若其他工具需要仿真时间，可显式传入 `use_sim_time:=true` 并用 `rosbag play --clock`；即使时钟暂停，节点的回调处理循环也不会因 `Rate::sleep()` 挂起。

按空格开始。首次建议半速，给球面插值和多相机事件处理留出计算余量；这是播放速度，不改变算法使用的消息时间戳。之后可按目标机吞吐能力改 `-r 1.0`，当前没有实时性测试结论。相邻分卷可按时间顺序传给同一个 `rosbag play`；重新播放一个序列时重启估计器，避免时间倒退。

相机数量可选：

```bash
# 私有数据：1 或 2 路；HILTI22 / M2DGR：1 至 5 路，按预设顺序取前 N 路。
roslaunch spherical_livo mapping_hilti22.launch num_cameras:=5
roslaunch spherical_livo mapping_m2dgr.launch num_cameras:=5
```

HILTI22 五路顺序为 **C0、C3、C4、C1、C2**，默认优先三个不同朝向；M2DGR 五路顺序为 **left、third、fourth、fifth、sixth**。M2DGR 不启用旧配置中外参尚属暂定的天空相机。`camera0…camera4` 是本项目连续编号，不等于原始硬件编号。

所有入口均使用同一共享视觉地图，并允许其他相机的历史参考参与当前相机更新，没有额外的跨相机开关。复用能否发生仍取决于共同可见区域、平面支持和有效光度观测。

### 私有数据的消息格式

默认采用旧配置对应的话题：`/livox/lidar`、`/livox/imu`、`/fisheye/left/image_raw`、`/fisheye/right/image_raw`。MEI 原图为 **1088 × 1280**，保留已有内外参，不先去畸变或缩放。

- 默认 `lidar_type: 8` 对应新增加的 Livox PointCloud2。上游编号 `7` 保留给 Robosense，不能照抄旧 cake_slam 的编号 7。
- PointCloud2 适配直接移植 cake_slam 的 `findField`、`readFieldAsDouble`、`normalizeLivoxPointTimeMs` 和 `livox_pointcloud2_handler`。时间字段按 `offset_time`、`timestamp`、`time`、`t` 的优先级选择；线号接受 `line/ring`，强度接受 `intensity/reflectivity`。时间减去有效点最小时间后，按整帧时间跨度使用旧项目的单位归一化规则。
- 若 bag 实际记录的是 `livox_ros_driver2/CustomMsg`，使用 `roslaunch spherical_livo mapping_private_mid360.launch lidar_type:=1`；原始 `livox_ros_driver/CustomMsg` 是不同 ROS 消息包，需先转换为 Driver 2 消息或带逐点时间的 PointCloud2。
- 缺少逐点时间，或逐点时间无有效跨度时，沿用 cake_slam 按点序和 `preprocess/scan_rate` 生成近似时间的处理，私有预设为 10 Hz；无时间字段时打印原有提示。这使输入兼容行为与旧项目一致，但生成的是近似时间，不能恢复已经丢失的真实采样时刻。ROS2 的 `.db3/.mcap` 仍需先转换成 ROS1 bag。
- 实机在线采集时沿用默认 `use_sim_time:=false`，由设备驱动发布上述话题；本项目启动文件只负责估计器。

### HILTI22 与 M2DGR 输入

HILTI22 使用 `/hesai/pandar` 的 PointCloud2（type 5，XT32）与 `/alphasense/imu`，图像来自 `/alphasense/cam*/image_raw`。已重新整理官方逐相机 KB4 内参；`Rcl/Pcl` 由官方 `T_cam_imu` 与 LiDAR→IMU 标定相乘得到，固定时间偏移使用各相机提供的 `timeshift_cam_imu`。若你使用的是已经修正图像时间戳的二次处理 bag，应把对应 `time_offset` 置零，避免重复补偿。[标定来源与编号](calibration_sources.md)

M2DGR 使用 `/velodyne_points` 与 `/handsfree/imu`；默认三路图像为原始 Image，额外两路为 `/camera/fifth/image_raw/compressed` 和 `/camera/sixth/image_raw/compressed`。压缩图像由节点内部解码，不需要启动 republish。Velodyne 沿用 FAST-LIVO2 预处理：有逐点 `time` 时按秒读取；缺失时使用上游按线束方位角、约 10 Hz 扫描恢复时间的近似，不能据此声称原始逐点时间已恢复。

如果你的 bag 话题命名或压缩形式不同，只需修改相应 `img_topic`/`lid_topic`/`imu_topic`。`/compressed` 后缀决定订阅的消息类型，后缀与实际 ROS 类型必须一致。自定义 YAML 可通过 `config:=/absolute/path.yaml cameras:=/absolute/intrinsics.yaml` 传入各数据集启动文件。

## 5. 默认参数为什么这样设

| 参数 | 私有 MID360 | HILTI22 | M2DGR |
| --- | --- | --- | --- |
| 球冠角半径 `patch_radius_deg` | 0.9° | 1.0° | 0.7° |
| 采样点数 `patch_samples` | 48 | 48 | 48 |
| 每次相机事件 patch 上限 | 120 | 120 | 120 |
| 点序抽稀 `point_filter_num` | 3 | 3 | 3 |
| 点云降采样 `filter_size_surf` | 0.10 m | 0.10 m | 0.10 m |
| 雷达盲区 | 0.80 m | 0.60 m | 0.60 m |
| 基础体素 `voxel_size` | 0.50 m | 0.50 m | 0.40 m |
| 平面最小特征值门限 | 0.0025 m² | 0.001 m² | 0.001 m² |
| IMU 初始化计数 | 100 | 100 | 100 |
| `acc_cov / gyr_cov` | 0.5 / 0.1 | 0.5 / 0.01 | 0.5 / 0.01 |

角半径依据中心角分辨率选取：MEI 的有效中心焦距约为 `f/(1+xi) ≈ 397 px/rad`，不是直接使用 1660；HILTI KB4 约为 `351 px/rad`；M2DGR 约为 `540 px/rad`。对应初始球冠半径约 6–7 个中心原始像素，48 个球面节点具有相近密度。这只是初始覆盖范围的比较，实际采样、warp 和带宽匹配仍在球面上，周边角分辨率使用 LUT 的局部信息。

点序抽稀从旧值 5 改为 3，保留更多雷达候选和平面／深度支持，再由 0.1 m 降采样与体素层级控制密度。HILTI22/M2DGR 的平面特征值门限从 0.0001 调至 0.001，减少对低噪声理想平面的过度要求；它对应法向散布标准差量级约 3.2 cm，并不取代深度一致性检查。私有 MID360 保留原有 0.0025。它们是静态工程起点，无法在不回放数据的情况下确认最优。

IMU 初始化计数在实现里按采样累计，并非相机帧数；提高到 100 为启动均值留出更多样本。应保留数据开头的静止段。IMU 噪声沿用代码的离散协方差约定，不能直接填入厂商的连续时间噪声密度。私有数据的 `gyr_cov` 从旧 0.8 改为较温和的 0.1，其余公开数据保留原配置量级；这不是新的 IMU 标定结果。

三级视觉尺度、插值支持、异常值处理和深度判据由实现统一处理，不需要恢复旧配置中的虚拟焦距、NCC 阈值、像素 patch、在线外参或方向性开关。开启的 `local_map/map_sliding_en` 是原有雷达地图空间裁剪，用于限制长序列地图范围，不是高级共享视觉地图管理。

## 6. 结果与常见部署问题

设置 `output_dir` 后生成：

- `trajectory.txt`：`timestamp tx ty tz qx qy qz qw`，输出 **IMU** 位姿。
- `visual.csv`：相机编号、有效 patch、跨相机参考计数、是否更新、残差及耗时。

ROS 输出为 `/spherical_livo/odometry`、`/spherical_livo/path`、`/spherical_livo/cloud`，固定坐标系 `world`。所有启动入口默认开启 RViz，共用 `rviz_cfg/spherical_livo.rviz`，已配置 `world` 坐标系、点云、轨迹和当前 IMU 位姿。点云显示累积最近 30 秒发布的数据，不是完整历史地图；显示累积不会影响算法地图。

已经部署的机器若尚未安装 RViz，执行 `sudo apt install ros-noetic-rviz`。无桌面环境时关闭可视化：

```bash
roslaunch spherical_livo mapping_private_mid360.launch rviz:=false
```

HILTI22、M2DGR 和通用 `mapping.launch` 使用相同开关。RViz 不是必需进程，关闭其窗口不停止估计器。源码工作空间拉取本次配置更新后无需重新编译 C++；使用 `catkin_make install` 的部署需重新安装，以复制新的配置文件。

| 现象 | 处理 |
| --- | --- |
| 找不到 `livox_ros_driver2` / `CustomMsg.h` | 先编译 ROS1 驱动，再 source 驱动工作空间后编译本包 |
| 编译进程 `cc1plus` 被 killed | 改为 `-j1 -l1`，通常是编译内存不足 |
| 找不到 `spherical_livo` 包 | 检查源码层级和 `source ~/spherical_ws/devel/setup.bash` |
| `Spherical LIVO ready` 前等待较久 | 正在建立原图 LUT，等准备完成再播放 bag |
| 图像尺寸不匹配 | 使用与内参相符的原图，不能把缩放或去畸变后的图像套用原始标定 |
| 提示 `has no time field` | 已启用旧适配层的按点序近似时间；`scan_rate` 是点云发布频率，私有预设为 10 Hz |
| 迟到图像或队列容量报错 | 先降低 bag 播放速度，并确认 IMU、点云、图像处于同一时间基准 |
| 有里程计但视觉 `updated` 长期为 0 | LiDAR 仍可独立运行；先查图像输入、标定和深度支持，不能当作视觉已正常融合 |

首次部署只需完成编译与一段数据回放，确认输入匹配及视觉实际更新。本文没有安排繁复检查，也没有把尚未执行的编译、精度评估或性能测试写成已通过。
