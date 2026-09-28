# Spherical LIVO — ROS1

在 FAST-LIVO2 原版主节点、IMU、LiDAR 和 VIOManager 上修改球面视觉观测。基准提交 `0d2c0346107b75b59934975adec9a6eeeb913c64`。

保留原版固定 19 维状态、LIO→VIO 顺序、体素地图、视觉地图生命周期、IEKF、着色发布和 ANSI 彩色计时表。多相机组帧和 MID360 PointCloud2 适配来自 cake_slam；多相机共享视觉地图并支持跨相机历史参考。新增球面几何不是另一个主节点或独立视觉估计器。

**本次仅做源码及基本静态核对，没有编译、运行测试、回放或提交。**

- [恢复来源及必要改动](docs/restoration_scope.md)
- [ROS1 编译部署](docs/deployment_ros1_zh.md)
- [球面模型及近似](docs/method.md)
- [标定来源](docs/calibration_sources.md)
- [研究 TODO](docs/research_todo.md)

```bash
source /opt/ros/noetic/setup.bash
source ~/livox_ros1_ws/devel/setup.bash
cd ~/spherical_ws
catkin_make -DCMAKE_BUILD_TYPE=Release -DCATKIN_ENABLE_TESTING=OFF -j2 -l2
source devel/setup.bash
roslaunch spherical_livo mapping_private_mid360.launch output_dir:=$HOME/spherical_results/run01
# 另一终端（已 source ROS 环境）
rosbag play /absolute/path/to/sequence.bag
```

普通回放不需要 `--clock`。默认打开共用 RViz，固定坐标系 `camera_init`，订阅 `/cloud_registered` 的 RGB8、`/path`、`/aft_mapped_to_init`、`/rgb_img/camera_0` 和 `camera_1`。其他相机可在 RViz 勾选。关闭显示用 `rviz:=false`。

私有 MID360、HILTI22、M2DGR 分别使用 `mapping_private_mid360.launch`、`mapping_hilti22.launch`、`mapping_m2dgr.launch`。保留既有 YAML 内外参；私有 PointCloud2 类型为 8，Livox CustomMsg 为 1。相机输入保持原始标定尺寸、scale=1，支持 MEI、KB4、Pinhole，以及原始／压缩图像。

恢复原版 Frame 后需要 Sophus、fmt、image_transport；已有 cake_slam/FAST-LIVO2 环境可复用，详见部署文档。源代码保留三项球面参数：`visual/patch_radius_deg`、`visual/patch_samples`、`visual/max_patches`。配置值尚未经回放调优。
