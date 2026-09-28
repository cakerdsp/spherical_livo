# 原版恢复范围

基准：本地 `upstream/main`，FAST-LIVO2 `0d2c0346107b75b59934975adec9a6eeeb913c64`。
恢复前的全部本地源码（包括未提交改动）另存为备份，未执行 git reset，未提交。

| 部件 | 当前来源与改动 |
| --- | --- |
| `main.cpp` | 恢复原版入口，仅改 ROS 节点名、私有参数命名空间 |
| `LIVMapper.h/.cpp` | 从原版恢复构造、初始化、IMU/LIO/VIO 顺序、handleLIO、地图更新、发布、保存、彩色计时表；组帧与多相机同步取自 cake_slam |
| `IMU_Processing.h/.cpp` | 恢复原版 19 维传播、偏置索引、积分与去畸变；仅补初始化消费区间的末端时间赋值、清空无效输出；time_list 定义移至 .cpp（同 cake_slam），补 OpenMP 头 |
| `voxel_map.h/.cpp` | 恢复该上游提交原文件；没有自写体素地图或 LIO 优化器 |
| `preprocess.h/.cpp` | 保留现有 FAST-LIVO2 预处理与 cake_slam MID360 PointCloud2 解析移植、Driver2 类型；没有撤销必需适配 |
| `common_lib.h` | 原版固定 19 维状态、加减与协方差；加多相机图像容器、类型 8、cake_slam 的 Sophus 兼容引用 |
| `vio.h/.cpp` | 从原版 VIOManager 修改。保留网格选点、VisualPoint/Feature 地图与观测生命周期、粗到细、7 列观测与 19 维 IEKF 更新公式、绘图和原版计时表。采样、warp、深度一致性和雅可比改为球面版本；不再调用独立 VisualEstimator |
| `frame.*`、`feature.h`、`visual_point.*` | 恢复原版类型、实现和所有权关系；增加源相机、球面图像／方向与有限平面支持字段 |
| `spherical.*`、`spherical_camera.h` | 必要新增的球面几何及原版 AbstractCamera 接口适配；不是独立前端或优化器 |

视觉检索还修正了原版负坐标处 floor 后重复减一与插入索引不一致的问题，保留插入侧约定。

删除独立 `spherical_vio.*`、`visualization.*`、`timing.h`，删除仅针对旧独立估计器的 integration test。保留球面几何测试源文件，但本次不运行。

主节点必要适配还包括：原配置键兼容、原始／压缩图像输入、各相机固定时间偏移、WallRate（普通 rosbag play 无需 --clock）、空点云保护、图像 LIO 已完成后立即进入 VIO、不重复等待已消费的雷达队列、切分点云时间排序、多相机图像／着色输出。既没有自写 IMU 端点插值，也没有另设一套异步调度。

多相机共享同一视觉地图和惯性状态，各相机顺序更新。逆曝光使用原版第 6 槽，切换相机时装入该相机上次曝光并重置该曝光的交叉协方差和先验方差；没有 20 维临时求解器。各图像组使用 cake_slam 的组时间（最大校正时间），容差内图像视为同一时刻；这不是连续时间异步多相机模型。

验证边界：仅基本源码、接口、配置与差异核对；没有编译、运行测试或回放。恢复原版框架不能替代部署验证。
