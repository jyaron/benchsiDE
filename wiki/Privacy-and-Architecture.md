# Privacy and architecture

## Where computation happens

benchsiDE is a single HTML file containing all application code and embedded resources (annotation tables, gene-set libraries, ortholog tables). When a file is loaded, it is read by the browser's File API into memory. All computation happens in the page. There is no server component and no upload.

## Network access

| Build | Requests |
|---|---|
| Hosted and `index.html` | The page itself and Plotly.js from its CDN |
| `benchside-offline.html` | None |

Gene links in the Gene Explorer open external databases in a new tab when clicked; they contain the gene identifier only.

## How to verify

Open the browser developer tools, select the Network tab, and run an analysis. No request carries expression data. With the offline build, the analysis can be run with networking disabled.

## Speed

Analyses complete in seconds because there is no upload, no queue and no server round trip, and the statistical methods (linear models with empirical-Bayes moderation) are closed-form. The speed does not come from approximating the statistics.
