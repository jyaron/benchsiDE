# Troubleshooting

| Symptom | Cause | Remedy |
|---|---|---|
| Every sample is in its own group | The design file was not applied. A common cause is a design saved with classic Mac line endings, which older versions read as one line. | Current versions read all line endings. Check that sample names in the design match the matrix headers exactly. |
| An extra "sample" named after the ID column | A numeric identifier column (Entrez IDs) was read as data | Current versions detect identifier columns. Place identifiers in the first column with a header such as `EntrezID`. |
| Genes appear as numbers or Ensembl IDs | No annotation applied | Choose the species and click **Use built-in annotation**, or load an annotation file. |
| A known gene is missing after annotation | The identifier is absent from the annotation release, or was filtered out | **Download unmapped IDs**; check filtering on the review step. |
| voom is not offered | The input is not raw counts | voom requires counts. Use the moderated t for normalized input. |
| Covariate option is absent | The design has one factor only | Add the covariate as a column in the design file. |
| "Covariate is confounded" message | The covariate separates the groups exactly | The effects cannot be separated. Redesign, or analyze without the covariate and state the confounding. |
| Few or no DE genes despite visible differences | Low replication; FDR unattainable | See [Statistical power](Statistical-Power.md). |
| Enrichment finds no sets | Wrong species library, or identifiers not decoded to symbols | Use the library matching the species; apply annotation first. |
| Cross-dataset match count is near zero | Different identifier types or species | Set the species of the added dataset; for other species load a **Custom ID map**. |
| Session loads with a version warning | The session was saved by another version | Results may differ where defaults changed; see CHANGELOG.md. |
| A plot looks crowded at poster size | Very long sample or group names | Labels are shortened automatically; shorter names in the design file give cleaner figures. |
| Self-test reports a failure | A numerical difference in this browser | Download the evidence file and open a GitHub issue with it. |
