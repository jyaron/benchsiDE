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

Every plot passes through a layout step that enlarges margins to fit titles, legends, annotations and tick labels at the current type size. Labels longer than the margin allows are shortened in the middle, keeping the start and the distinguishing end (for example `Mouse_skin…PMR_1_count`); the full name remains in the hover text. If margins would reduce the plotting area below 160 pixels, the figure grows instead.

## Composite figures

| Export | Location | Content |
|---|---|---|
| Summary figure | Header | Four-panel overview with caption |
| Panel figure | Gene Explorer, Discovery, Enrichment | One panel per gene with pairwise statistics; caption file |
| Module heatmap | Discovery, Enrichment | All member genes with a module-score strip |

Exported figures contain no interactive elements or on-screen instructions.

## Tables

Most panels offer CSV export of the full underlying table, not only the rows displayed.
