# Privacy and architecture

## Where computation happens

benchsiDE is a single HTML file containing all application code and embedded resources (annotation tables, gene-set libraries, ortholog tables). When a file is loaded, it is read by the browser's File API into memory. All computation happens in the page. There is no server component and no upload.

## Network access

| Build | Requests |
|---|---|
| Hosted and `index.html` | The page itself and Plotly.js from its CDN |
| `benchside-offline.html` | None |

Gene links in the Gene Explorer open external databases in a new tab when clicked; they contain the gene identifier only.

## Content-Security-Policy

Both builds carry a Content-Security-Policy, enforced by the browser, that forbids the page from opening network connections (`connect-src 'none'`: no fetch, XMLHttpRequest, WebSocket or beacon) and from loading images, scripts, styles, fonts, frames or form targets from any server. The single exception is Plotly.js in the online build, loaded from `https://cdn.plot.ly` and pinned by a subresource-integrity hash, so a modified copy would not run. The offline build allows no server at all. The policy is in the page source (`<meta http-equiv="Content-Security-Policy" ...>`) and cannot be relaxed by the page's own code.

Two allowances are required by the software and do not permit network access: `'unsafe-inline'`, because the application is one inline script, and `'unsafe-eval'`, because Plotly's WebGL scatter plots compile their shaders with the Function constructor.

The policy does not govern links the user clicks (gene pages, licence pages), which open in a new tab, or files the user saves.

## How to verify

1. Open the browser developer tools and the Network tab, then run an analysis. The only request is the page itself (and Plotly.js for the online build).
2. In the Console, type `fetch("https://example.org")`. The browser refuses it and reports a Content-Security-Policy violation.
3. With the offline build, the analysis can also be run with networking switched off.

The continuous-integration test `tests/privacy.spec.js` performs these checks in Chromium, Firefox and WebKit on both builds.

## Speed

Analyses complete in seconds because there is no upload, no queue and no server round trip, and the statistical methods (linear models with empirical-Bayes moderation) are closed-form. The speed does not come from approximating the statistics.
