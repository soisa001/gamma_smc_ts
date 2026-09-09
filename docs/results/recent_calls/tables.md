# Audited tables: posterior calls and selection onset

Every rate cell is **selected regions called / neutral regions called** (100 selected and 1,000 shared neutral regions). A call anywhere in the isolated 10 Mb region counts. All selected replicates are retained.

## Allele frequencies

| Onset (kya) | n | Mean | Median | Min | Max | Fixed | AF < 0.1 |
|---|---|---|---|---|---|---|---|
| 50 | 100 | 0.783925 | 0.82125 | 0.2575 | 0.9925 | 0 | 0 |
| 10 | 100 | 0.212400 | 0.17250 | 0.0075 | 0.6600 | 0 | 31 |

## Regional maximum at nominal alpha 0.05

Selection onset **50 kya**.

| T (kya) | mean | median | prob80 | prob90 | truth |
|---|---|---|---|---|---|
| 5 | 12% / 4.2% | 17% / 4.7% | 3% / 3.5% | 3% / 4.1% | 1% / 3.4% |
| 10 | 18% / 4.5% | 18% / 5.8% | 18% / 4.6% | 7% / 4.2% | 2% / 4.3% |
| 20 | 19% / 5.5% | 16% / 5.3% | 17% / 5.2% | 18% / 4.6% | 9% / 5.4% |
| 30 | 19% / 5.5% | 13% / 5.2% | 18% / 5.1% | 16% / 5.4% | 10% / 4.5% |
| 40 | 16% / 5.4% | 12% / 4.5% | 18% / 5.3% | 18% / 4.6% | 24% / 4.9% |
| 50 | 14% / 4.9% | 10% / 4.2% | 17% / 5.3% | 19% / 5.4% | 45% / 5.5% |

Selection onset **10 kya**.

| T (kya) | mean | median | prob80 | prob90 | truth |
|---|---|---|---|---|---|
| 5 | 7% / 4.2% | 3% / 4.7% | 2% / 3.5% | 3% / 4.1% | 1% / 3.4% |
| 10 | 3% / 4.5% | 3% / 5.8% | 6% / 4.6% | 3% / 4.2% | 8% / 4.3% |
| 20 | 5% / 5.5% | 7% / 5.3% | 2% / 5.2% | 2% / 4.6% | 5% / 5.4% |
| 30 | 8% / 5.5% | 7% / 5.2% | 5% / 5.1% | 3% / 5.4% | 6% / 4.5% |
| 40 | 7% / 5.4% | 4% / 4.5% | 4% / 5.3% | 5% / 4.6% | 6% / 4.9% |
| 50 | 6% / 4.9% | 5% / 4.2% | 7% / 5.3% | 6% / 5.4% | 5% / 5.5% |

## Regional maximum at nominal alpha 0.1

Selection onset **50 kya**.

| T (kya) | mean | median | prob80 | prob90 | truth |
|---|---|---|---|---|---|
| 5 | 17% / 8.7% | 25% / 10.0% | 5% / 7.7% | 6% / 9.2% | 5% / 9.1% |
| 10 | 24% / 10.2% | 22% / 9.9% | 29% / 10.2% | 10% / 9.3% | 5% / 8.8% |
| 20 | 28% / 10.9% | 26% / 9.7% | 24% / 9.7% | 28% / 10.2% | 19% / 10.4% |
| 30 | 27% / 10.2% | 25% / 9.7% | 27% / 11.3% | 25% / 10.1% | 18% / 8.9% |
| 40 | 24% / 10.3% | 21% / 10.2% | 28% / 9.8% | 27% / 11.1% | 31% / 9.6% |
| 50 | 23% / 9.8% | 18% / 9.5% | 28% / 10.3% | 26% / 11.1% | 47% / 10.6% |

Selection onset **10 kya**.

| T (kya) | mean | median | prob80 | prob90 | truth |
|---|---|---|---|---|---|
| 5 | 8% / 8.7% | 13% / 10.0% | 5% / 7.7% | 5% / 9.2% | 7% / 9.1% |
| 10 | 8% / 10.2% | 6% / 9.9% | 10% / 10.2% | 8% / 9.3% | 10% / 8.8% |
| 20 | 9% / 10.9% | 13% / 9.7% | 9% / 9.7% | 11% / 10.2% | 9% / 10.4% |
| 30 | 12% / 10.2% | 10% / 9.7% | 10% / 11.3% | 9% / 10.1% | 12% / 8.9% |
| 40 | 10% / 10.3% | 11% / 10.2% | 12% / 9.8% | 13% / 11.1% | 12% / 9.6% |
| 50 | 10% / 9.8% | 12% / 9.5% | 11% / 10.3% | 10% / 11.1% | 10% / 10.6% |

## Thresholds trained for 70% selected sensitivity

These fixed T=10 and T=50 kya examples use the raw regional maximum. Each threshold is trained on 80 selected regions, then evaluated on held-out regions. The table reports realized held-out rates, which need not equal the training target. All other cutoffs and the standardized score remain available in region_calls.csv.

| Rule | T (kya) | Onset 50 kya | Onset 10 kya |
|---|---|---|---|
| mean | 10 | 68% / 54.0% | 68% / 67.9% |
| mean | 50 | 71% / 57.1% | 69% / 68.5% |
| median | 10 | 69% / 56.5% | 71% / 69.2% |
| median | 50 | 71% / 59.6% | 67% / 69.0% |
| prob80 | 10 | 69% / 48.8% | 68% / 65.4% |
| prob80 | 50 | 70% / 54.7% | 68% / 66.8% |
| prob90 | 10 | 73% / 65.8% | 70% / 73.9% |
| prob90 | 50 | 68% / 59.8% | 70% / 64.6% |
| truth | 10 | 73% / 70.5% | 77% / 66.3% |
| truth | 50 | 68% / 25.9% | 70% / 69.9% |

## Stride p <= 0.001 with consecutive evidence

These examples use 10 kb reporting strides, with five or ten significant strides required and no gaps. The pointwise p threshold does not imply a regional 0.001 error rate. The full exploratory grid is in paired_grid.csv; no winner from that grid has received separate confirmation.

TMRCA cutoff **10 kya**.

| Rule | Required strides | Onset 50 kya | Onset 10 kya |
|---|---|---|---|
| mean | 5 | 25% / 7.7% | 7% / 7.7% |
| mean | 10 | 17% / 3.0% | 1% / 3.0% |
| median | 5 | 24% / 7.4% | 7% / 7.4% |
| median | 10 | 12% / 1.6% | 1% / 1.6% |
| prob80 | 5 | 20% / 6.5% | 6% / 6.5% |
| prob80 | 10 | 12% / 2.8% | 2% / 2.8% |
| prob90 | 5 | 4% / 4.5% | 4% / 4.5% |
| prob90 | 10 | 3% / 2.7% | 4% / 2.7% |
| truth | 5 | 1% / 4.0% | 5% / 4.0% |
| truth | 10 | 0% / 1.0% | 2% / 1.0% |

TMRCA cutoff **50 kya**.

| Rule | Required strides | Onset 50 kya | Onset 10 kya |
|---|---|---|---|
| mean | 5 | 17% / 3.4% | 1% / 3.4% |
| mean | 10 | 3% / 0.3% | 0% / 0.3% |
| median | 5 | 7% / 1.8% | 2% / 1.8% |
| median | 10 | 1% / 0.1% | 0% / 0.1% |
| prob80 | 5 | 22% / 5.7% | 0% / 5.7% |
| prob80 | 10 | 8% / 0.7% | 0% / 0.7% |
| prob90 | 5 | 22% / 7.6% | 6% / 7.6% |
| prob90 | 10 | 10% / 1.0% | 0% / 1.0% |
| truth | 5 | 34% / 3.7% | 4% / 3.7% |
| truth | 10 | 15% / 0.3% | 0% / 0.3% |

## Provenance

Analysis source: `D:\phase2simselection\sim\eas_recent_calls\analysis`. Input/output inventory SHA-256: `06aa76baddffe53a3a3ab37bdb64758351b603d613406beede56549fbad9f922`. Table-rendering script SHA-256: `63f6baa10ef650c62745e227f703c5ec284d09474cf0a1d1ed8c1d8e8a98ce4a`. The independent audit recomputed all 960 regional summaries from 528,000 unique held-out predictions and reproduced the 192 existing mean/truth baseline rows.
