# Orthologs and homology families

## One-to-one orthologs

Mouse–human gene matching uses the MGI homology report (HOM_MouseHumanSequence.rpt), restricted to strict one-to-one pairs. This resolves renamed orthologs such as Trp53 and TP53, which symbol matching misses.

## Genes without one-to-one orthologs

Some gene families expanded independently in each lineage. The mouse genes Serpinb3a, Serpinb3b, Serpinb3c and Serpinb3d correspond jointly to human SERPINB3 and SERPINB4, with no one-to-one pairing. Such genes are excluded from gene-level matching, because any pairing would be arbitrary.

## Homology families

To retain these genes, benchsiDE builds homology families: connected groups of the homology classes, merged by shared members (union–find), so that a whole expanded family forms one unit.

**Score families** (Compare datasets tab) treats each family as a module:

1. In each dataset separately, the family score is the mean z-score of its members present in that dataset.
2. The group effect is tested within each dataset by Welch t on the family scores, with BH correction across families.
3. Agreement between datasets is reported as direction concordance.

## Family concordance scatter

One point per family: x is the effect in the session, y the effect in the comparison dataset. Point size reflects the number of matched members.

| Colour | Meaning |
|---|---|
| Red | Concordant up; significant in both datasets |
| Blue | Concordant down; significant in both datasets |
| Amber diamond | Significant in both datasets, opposite directions |
| Grey | Not significant in at least one dataset |

Click a point, or **plot** in the table, for per-sample scores in both datasets.

## Custom mappings

For other species pairs, or to override the defaults, load a two-column **Custom ID map** (source identifier, target identifier). Identifiers mapping to more than one target are dropped.
