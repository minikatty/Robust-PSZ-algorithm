# Robust-PSZ-algorithm
The implementation of a robust hybrid method for PSZ

## 1 Project structure

```
Project/                         %项目主目录
│
├── README.md                    % 项目概览说明文档
├── main.m                       % 主程序入口
├── config/                      % 参数配置文件夹
│   ├── config_default.mat       % 默认参数（fs, beta, roomSize 等）
│   └── config_custom.mat        % 具体实验配置
│
├── RIRdata/                     % 输入 / 输出数据目录
│   ├── raw/                     % 原始数据（测量/录制）
│   ├── processed/               % 处理后结果（滤波、降采样）
│   └── results/                 % 实验输出结果（RIR、IR等）
│
├── src/                         % 核心源代码
│   ├── rir_generator.m          % RIR生成函数
│   ├── rir_postprocess.m        % RIR预处理（滤波、截断、归一化）
│   ├── showStruct.m             % 自动结构体分析工具（见下）
│   └── utils/                   % 工具函数（如绘图、分析等）
│
├── docs/                        % 文档、图示、流程图等
│   ├── structure_diagram.png    % 数据结构图
│   └── signal_flow_mermaid.md   % Markdown信号流图
│
└── logs/                        % 日志、运行记录
    └── meta_2025_10_04.mat      % 运行时参数记录

```


##  2 data structure

```
IR
 ├─ HB_ctrl      (cell) 控制点RIR
 ├─ HB_eval      (cell) 评价点RIR
 ├─ parameters   (struct)
 │   ├─ fs       (double) 采样率
 │   ├─ beta     (double) 混响时间
 │   └─ roomSize (1x3 double) 房间尺寸

```

