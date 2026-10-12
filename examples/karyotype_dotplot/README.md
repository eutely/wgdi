# Karyotype dotplot：核型共线性点图（`-kd`）

本模块采用 WGDI 原生输入和 Matplotlib 绘图方式。色块直接来自核型文件，
黑点独立来自共线性结果；色块的范围和颜色不随共线性点的密度变化。
无需新增安装依赖。

```bash
# 在本示例目录运行
wgdi -kd example.conf
# 导出完整配置模板
wgdi -kd '?' > karyotype_dotplot.conf
```

## 输入方向：与 WGDI dotplot 一致

| 参数 | 含义 |
| --- | --- |
| `gff1`, `lens1` | 现存物种（left），显示为竖直染色体，lens 行顺序决定排列 |
| `gff2`, `lens2` | 祖先核型（top），决定黑点在每条染色体内部的横坐标 |
| `collinearity` | WGDI `-icl` 输出，第一组基因对应 gff1，第二组对应 gff2 |
| `ancestor_left` | 现存物种的核型分区文件，决定全部背景色块 |
| `ancestor_top` | 祖先的核型文件，用于将图例颜色标成祖先染色体编号；可省略或设为 none |
| `genome1_name`, `genome2_name` | 现存物种轴标签与祖先图例标题 |

`ancestor_left` 也可使用 WGDI `-k` 的参数名 `ancestor`；两者若同时提供，
必须指向同一个文件。所有相对文件路径以运行命令的工作目录为基准。

### GFF 与 lens

使用 WGDI 的简化 GFF，**不是 GFF3**。制表符分隔，可保留额外列：

```text
chromosome  gene_id  start  end  strand  order
```

其中 order 是染色体内从 1 开始的基因序号。每个基因 ID 唯一，每条染色体内
order 唯一。lens 为三列：

```text
chromosome  bp_length  gene_count
```

lens 同时选择染色体、指定显示顺序和完整长度。程序复用 `base.newgff()`、
`base.newlens()`，并检查长度是否覆盖注释坐标。

### 共线性文件

直接使用 `wgdi -icl` 的输出，不需要转换或交换输入列：

```text
# Alignment 1: score=100 pvalue=0.01 N=3 s1&a1 plus
species_gene1 1 ancestor_gene1 1 1
species_gene2 2 ancestor_gene2 2 1
species_gene3 3 ancestor_gene3 3 1
```

复用 `base.read_collinearity()`。默认严格检查基因 ID 与输入方向；不会自动交换
列。没有匹配到的基因会报错，设 `strict = false` 可改为警告后跳过。

### 核型文件：直接指定颜色及完整范围

与 WGDI `-k`、`-km`、`ancestor_left/ancestor_top` 相同，制表符分隔五列：

```text
chromosome  start_order  end_order  color  weight
```

例如：

```text
1	1	807	RoyalBlue	1
1	808	2266	fuchsia	1
2	1	1474	#FFCC00	1
```

- 第 2、3 列是 **从 1 开始、两端包含的基因序号**，不是 bp。
- 第 4 列直接作为 Matplotlib 颜色，不重新分配调色板。
- 第 5 列沿用 WGDI 绘图行为，除以所选核型区间的最大值后作为透明度。
- 使用 `base.read_classification()` 和 `base.Rectangle()` 读取并绘制。
- 对每条 lens 选中的染色体，核型区间必须完整覆盖 `1..gene_count`，不能有
  缺口或重叠。未提供区间的染色体也会报错，避免静默编造祖先归属。
- 即使某区间没有共线性基因对，也会完整着色。程序不根据共线块做插值、
  平滑或补色。透明颜色或第 5 列的零值仍按输入绘制。

在 `position = order` 下，核型区间 `[start, end]` 绘制为 `[start-1, end]`，
因此相邻区间无缝衔接，首尾覆盖整个染色体。

在 `position = end` 下，核型文件仍按基因序号读取。内部边界放在前一区间末基因
与后一区间首基因的 **end 坐标中点**；首尾延伸到 0 和 lens 中的物理长度。
相邻区间共享同一边界，因此物理坐标图也连续覆盖所有位置。此换算表示绘图分界，
不是精确推断的重排断点。缺少边界基因或换算后的边界不递增会报错。

## 黑点坐标

沿用 WGDI 的 `position` 含义：`order` 使用原始一基基因序号，`end` 使用 GFF
的基因终点，不使用区间中点。设现存染色体的图中中心位置为 i：

```text
x = i - width/2 + (ancestor_gene_position / ancestor_chromosome_length) * width
y = extant_gene_position / maximum_extant_chromosome_length
```

所以点阵确实比较祖先与现存物种。改变核型着色文件只影响背景色块，
不改变基因对或黑点坐标。

## 配置示例

```ini
[karyotype_dotplot]
gff1 = species.gff
gff2 = ancestor.gff
lens1 = species.lens
lens2 = ancestor.lens
collinearity = species.ancestor.collinearity
ancestor_left = species.ancestor.txt
ancestor_top = ancestor.ancestor.txt
genome1_name = Extant species
genome2_name = Ancestral karyotype
position = order
width = 0.6
markersize = 0.5
figsize = 14,6.18
savefig = species.ancestor.svg
```

| 可选参数 | 默认值 | 说明 |
| --- | --- | --- |
| `position` | order | order / end，与 WGDI 相同 |
| `width` | 0.6 | 染色体宽度，相邻中心距离为 1；范围 (0,1] |
| `markersize` | 0.5 | Matplotlib 散点面积，单位 pt² |
| `markeralpha` | 0.5 | 黑点透明度，范围 (0,1] |
| `figsize` | 10,6.18 | 图宽、高，单位英寸 |
| `fontsize` | 10 | 标签字号 |
| `dpi` | 500 | 输出分辨率 |
| `block_length` | 1 | 共线块最小基因对数，只筛选黑点，不改变色块 |
| `labels`, `legend` | true | 显示染色体标签、祖先颜色图例 |
| `strict` | true | 缺失或方向不匹配的基因是否报错 |
| `title` | 空 | 图标题 |
| `savefig` | wgdi.karyotype_dotplot.svg | 支持 .svg、.pdf、.png |

## 与早期移植版本的区别

模块现命名为 `karyotype_dotplot`，命令为 `-kd`，配置节为
`[karyotype_dotplot]`。旧名称 `karyotype_svg` / `-ksvg` 已替换，
请同步修改命令及配置节名。

当前模块专门采用 WGDI 原生输入，已移除 JCVI BED/anchors、`mode`、`swap`、
`auto_swap`、`color_file`、`plot_style`、`band_style`、`max_gap` 等早期接口。
需要 JCVI 绘图时仍可使用工作区中原有的 Perl 脚本。

迁移配置时必须核对 **gff1 现存物种 / gff2 祖先** 的方向，使用 `position`
替换 mode，使用 `ancestor_left` 提供完整现存染色体分区。
原先 `ancestor` 指向 query 颜色区间的含义已取消；现在它仅是
`ancestor_left` 的别名，与 WGDI `-k` 一致。

## 验证

从仓库根目录运行：

```bash
MPLBACKEND=Agg python -m unittest discover -s tests -v
```

测试覆盖原生文件输入、完整区间覆盖、无点区间着色、透明度、坐标与方向、
端点换算、输入错误、SVG/PDF/PNG 及命令入口。
