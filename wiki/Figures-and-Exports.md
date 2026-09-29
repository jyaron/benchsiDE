# Figures and exports

## Plot export

Each plot has a toolbar (visible on hover). The camera icon saves the plot as SVG. The **export** selector sets the aspect ratio of exported figures: as shown, single column (4:3 or 1:1), or double column (7:4 or wide).

## Type size

| Setting | Scale | Use |
|---|---|---|
| Screen fonts | 1× | On-screen analysis |
| Paper fonts | 1.5× | Journal figures |
| Poster fonts | 2× | Posters and slides |

All plots, including exports, use the chosen scale.

## Layout

Every plot passes through a layout step that enlarges margins to fit titles, legends, annotations and tick labels at the current type size. If margins would reduce the plotting area below 160 pixels, the figure grows instead.

Axis tick labels are controlled from the header:

| Control | Options |
|---|---|
| Orientation | auto, horizontal, 45°, vertical. Auto uses horizontal labels when they fit, wraps them onto up to three lines, and otherwise angles them. |
| Label size | auto, or a fixed size from 8 to 16 (scaled by the type-size setting) |
| Long labels | wrap (default), full length, or shortened. Wrap and full length never cut a label off; the margin grows instead. Shortened labels keep the start and the distinguishing end (for example `Mouse_skin…PMR_1_count`), and the full name remains in the hover text. |

## Data for re-plotting (Excel / GraphPad Prism)

The **Download data** button in each plot's toolbar saves an Excel workbook of the plotted values:

| Sheet | Content |
|---|---|
| Prism | The data in GraphPad Prism table layout: Column (one column per group, replicate values stacked), Grouped (row titles plus one column per series; heatmaps as a matrix; forest plots with lower and upper limits) or XY (X column and one Y column per series). |
| Sample IDs | For column tables, the sample behind each value, in the same positions |
| All values | Every plotted point with series, label, coordinates and error values |
| Notes | Figure and axis titles, the Prism table type, import instructions, application version and export time |

In Prism, create a table of the stated type and paste or import the Prism sheet with the first row as column titles. The Prism sheet contains individual values rather than means, so that Prism computes its own summaries.

## Composite figures

| Export | Location | Content |
|---|---|---|
| Summary figure | Header | Four-panel overview with caption |
| Panel figure | Gene Explorer, Discovery, Enrichment | One panel per gene with pairwise statistics; caption file |
| Module heatmap | Discovery, Enrichment | All member genes with a module-score strip |

Exported figures contain no interactive elements or on-screen instructions.

## Tables

Most panels offer CSV export of the full underlying table, not only the rows displayed.
