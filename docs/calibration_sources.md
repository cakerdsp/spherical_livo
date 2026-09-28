# 数据集标定来源

配置整理日期：2026-09-27。算法参数和实测标定分开处理：前者给出未回放的默认值，后者保留数值来源，不以“调参”方式改造外参。

## 私有 MID360 + MEI

来自本地 cake_slam 的活动配置：

- `config/mid360_mei_multi_cam.yaml`：两个相机的 `Rcl/Pcl`、LiDAR→IMU 外参和输入话题。
- `config/camera_mei_multi_cam.yaml`：原始 camera0/1，1088 × 1280，MEI 模型。忽略注释的旧标定与 `camera0_undistort`。

固定时间偏移保留为零，因为没有提供其他已标定数值。未重新标定，未核实具体 bag 是否使用了相同传感器安装与图像格式。

## HILTI22

原 cake_slam 多相机内参文件明确注明把单个相机的内参复制到所有视角。本项目改为官方逐相机标定，来源是 [Hilti-Research 官方校准目录](https://huggingface.co/datasets/Hilti-Research/hilti-slam-challenge-2022/tree/main/calibration/calibration_files)，上传记录 `6ed5fe2`：

| 本项目编号 | 硬件编号 | 官方文件 | 加到图像时间戳的秒数 |
| --- | --- | --- | --- |
| camera0 | C0 | [calib_3_cam0-1-camchain-imucam.yaml](https://huggingface.co/datasets/Hilti-Research/hilti-slam-challenge-2022/blob/main/calibration/calibration_files/calib_3_cam0-1-camchain-imucam.yaml)，cam0 节 | 0.0019071204237431304 |
| camera1 | C3 | [calib_3_cam3-camchain-imucam.yaml](https://huggingface.co/datasets/Hilti-Research/hilti-slam-challenge-2022/blob/main/calibration/calibration_files/calib_3_cam3-camchain-imucam.yaml) | 0.0018128742099424312 |
| camera2 | C4 | [calib_3_cam4-camchain-imucam.yaml](https://huggingface.co/datasets/Hilti-Research/hilti-slam-challenge-2022/blob/main/calibration/calibration_files/calib_3_cam4-camchain-imucam.yaml) | 0.0017123033924705262 |
| camera3 | C1 | 上述 cam0-1 文件，cam1 节 | 0.001901867458041627 |
| camera4 | C2 | [calib_3_cam2-camchain-imucam.yaml](https://huggingface.co/datasets/Hilti-Research/hilti-slam-challenge-2022/blob/main/calibration/calibration_files/calib_3_cam2-camchain-imucam.yaml) | 0.001939579925335353 |

官方的 `pinhole + equidistant` 对应本实现的 `EquidistantCamera`（KB4），不能加载成普通无畸变针孔。单相机校准文件的顶层名字也可能叫 `cam0`，实际传感器身份按其中的 `rostopic` 确定。

LiDAR→IMU 由 [lidar_calibration.yaml](https://huggingface.co/datasets/Hilti-Research/hilti-slam-challenge-2022/blob/main/calibration/calibration_files/lidar_calibration.yaml) 给出：平移 `[-0.001, -0.00855, 0.055] m`；四元数为 xyzw `[0.7071068, -0.7071068, 0, 0]`，归一化后对应 `[[0,-1,0],[-1,0,0],[0,0,-1]]`。

采用 `p_C = R_CI p_I + t_CI` 和 `p_I = R_IL p_L + t_IL`，生成 `Rcl = R_CI R_IL`、`Pcl = R_CI t_IL + t_CI`。时间约定按 [Kalibr 文档](https://github.com/ethz-asl/kalibr/wiki/yaml-formats)：`t_imu = t_cam + timeshift_cam_imu`，无需在线标定。

校准数值的原作者为 Hilti-Research / Hilti-Oxford Dataset 团队，来源目录声明 **CC BY-NC-SA 3.0**；本项目将其重排并转换外参表示，未重新估计标定。`config/hilti22.yaml` 与 `config/cameras_hilti22.yaml` 中派生的标定内容遵循该许可：[许可全文](https://creativecommons.org/licenses/by-nc-sa/3.0/)。算法源码继续使用根目录 GPL-2.0。

## M2DGR

沿用此前从 cake_slam 的 `config/M2DGR_multi_cam.yaml` 与 `config/camera_M2DGR_multi_cam.yaml` 迁入的内外参和话题。五路顺序为 left、third、fourth、fifth、sixth；不包含暂定天空相机。未声称这些标定已在本实现回放验证。配置里的噪声和几何门限是默认参数，不是额外实测标定。
