# DEXSeq's DTU Model: How the Interaction Term Tests Usage

This note explains how the **DEXSeq** differential usage model works, with special focus on the **interaction term**.

It also explains how this differs from the simpler offset model described in [glm_dtu_offset_explanation.md](/home/unix/bhaas/projects/GITHUB/brianjohnhaas/GLMs_explained-/glm_dtu_offset_explanation.md), which is closer to the logic used by **saseR**.

The short version is:

> **saseR-style offset model:** put the gene total into the model as a fixed offset.  
> **DEXSeq interaction model:** put both the feature count and the other-gene count into the model, then use a condition-by-feature interaction to test whether their ratio changes.

Both are trying to answer the same biological question:

> Did this feature's usage change relative to the rest of its gene?

But they get there through different math.

---

## 1. The same toy gene

Suppose we have one gene with two transcript isoforms or exon bins:

- Feature **F1**
- Feature **F2**

We compare two groups:

- **Control**
- **Case**

Here are six samples:

| sample | condition | gene total counts | F1 counts | F2 counts |
|---|---|---:|---:|---:|
| C1 | control | 100 | 80 | 20 |
| C2 | control | 120 | 96 | 24 |
| C3 | control | 80 | 64 | 16 |
| D1 | case | 200 | 80 | 120 |
| D2 | case | 180 | 72 | 108 |
| D3 | case | 220 | 88 | 132 |

The gene total is:

$$
G_i = Y_{F1,i} + Y_{F2,i}
$$

For sample C1:

$$
G_{C1} = 80 + 20 = 100
$$

For sample D1:

$$
G_{D1} = 80 + 120 = 200
$$

Average usage values are:

| feature | control usage | case usage |
|---|---:|---:|
| F1 | 80% | 40% |
| F2 | 20% | 60% |

So the gene changed in two ways:

1. The total gene expression increased from 100 to 200.
2. The feature usage shifted from mostly F1 to mostly F2.

The DTU question is about the second point.

---

## 2. First recall the saseR-style offset model

The offset tutorial describes a model like this:

$$
Y_{f i} \sim NB(\mu_{f i}, \phi_f)
$$

$$
\log(\mu_{f i}) = \log(G_i) + \alpha_f + \beta_f x_i
$$

where:

| term | meaning |
|---|---|
| $Y_{f i}$ | observed count for feature $f$ in sample $i$ |
| $\mu_{f i}$ | expected count for feature $f$ in sample $i$ |
| $G_i$ | total parent-gene count in sample $i$ |
| $\log(G_i)$ | gene-total offset |
| $x_i$ | condition indicator: 0 for control, 1 for case |
| $\alpha_f$ | baseline log usage of feature $f$ |
| $\beta_f$ | condition effect on feature usage |

The important algebra is:

$$
\log(\mu_{f i}) = \log(G_i) + \alpha_f + \beta_f x_i
$$

Subtract the offset:

$$
\log(\mu_{f i}) - \log(G_i) = \alpha_f + \beta_f x_i
$$

Use the log ratio rule:

$$
\log\left(\frac{\mu_{f i}}{G_i}\right) = \alpha_f + \beta_f x_i
$$

So the model is really saying:

> log expected feature usage = baseline usage + condition effect

For F1:

$$
\beta_{F1,\ offset} = \log\left(\frac{0.4}{0.8}\right) = \log(0.5) \approx -0.693
$$

For F2:

$$
\beta_{F2,\ offset} = \log\left(\frac{0.6}{0.2}\right) = \log(3) \approx 1.099
$$

In this model, $\beta_f$ is a **log usage ratio**.

For F1:

$$
e^{\beta_{F1,\ offset}} = 0.5
$$

That means:

> F1 usage in case is 0.5 times F1 usage in control.

For F2:

$$
e^{\beta_{F2,\ offset}} = 3
$$

That means:

> F2 usage in case is 3 times F2 usage in control.

This is the model most relevant to saseR's adapted-offset idea.

---

## 3. DEXSeq uses a different setup

DEXSeq does not start by putting $\log(G_i)$ into the model as an offset.

Instead, for each feature being tested, DEXSeq compares:

1. the count for **this feature**, and
2. the count for **the rest of the gene**.

For F1, DEXSeq builds a table like this:

| sample | condition | row type | count |
|---|---|---|---:|
| C1 | control | this feature, F1 | 80 |
| C1 | control | other features | 20 |
| C2 | control | this feature, F1 | 96 |
| C2 | control | other features | 24 |
| C3 | control | this feature, F1 | 64 |
| C3 | control | other features | 16 |
| D1 | case | this feature, F1 | 80 |
| D1 | case | other features | 120 |
| D2 | case | this feature, F1 | 72 |
| D2 | case | other features | 108 |
| D3 | case | this feature, F1 | 88 |
| D3 | case | other features | 132 |

For each original RNA-seq sample, there are now two model rows:

- one row for the feature under test,
- one row for everything else in the gene.

This lets DEXSeq ask:

> Did the ratio of "this feature" to "the rest of the gene" change between conditions?

That ratio is closely related to usage.

For F1 in controls:

$$
\frac{\text{F1}}{\text{other}} = \frac{80}{20} = 4
$$

For F1 in cases:

$$
\frac{\text{F1}}{\text{other}} = \frac{80}{120} = \frac{2}{3} \approx 0.667
$$

So F1 went from 4 times the rest of the gene to only two-thirds of the rest of the gene.

That is a strong usage change.

---

## 4. The DEXSeq formula

In the DEXSeq R interface, the model is often written as:

```r
~ sample + exon + condition:exon
```

The reduced model is:

```r
~ sample + exon
```

The extra term in the full model is:

```r
condition:exon
```

That is the interaction term.

It means:

> Does the effect of being the feature row versus the other-gene row depend on condition?

That sentence is dense, so we will unpack it slowly.

---

## 5. What is an interaction?

An interaction means:

> The effect of one variable changes depending on another variable.

Here the two variables are:

| variable | values | question it asks |
|---|---|---|
| condition | control or case | Which biological group is this sample in? |
| exon row | this feature or other features | Are we counting this feature or the rest of the gene? |

Without an interaction, the model can say:

> Cases have more counts overall than controls.

and:

> This feature has more or fewer counts than the rest of the gene.

But without the interaction, it cannot say:

> This feature changed differently from the rest of the gene in cases.

That last sentence is exactly the DTU signal.

---

## 6. Define a tiny DEXSeq model by hand

For one feature, define:

$$
r = \text{row type}
$$

There are two row types:

$$
r =
\begin{cases}
0 & \text{other features} \\
1 & \text{this feature}
\end{cases}
$$

Define condition:

$$
x_i =
\begin{cases}
0 & \text{control} \\
1 & \text{case}
\end{cases}
$$

Let:

$$
Y_{i r} = \text{observed count for sample } i \text{ and row type } r
$$

Let:

$$
\mu_{i r} = \text{expected count for sample } i \text{ and row type } r
$$

A simplified DEXSeq-style model is:

$$
Y_{i r} \sim NB(\mu_{i r}, \phi)
$$

$$
\log(\mu_{i r}) = s_i + \alpha r + \beta x_i r
$$

Term by term:

| term | meaning |
|---|---|
| $s_i$ | sample effect for sample $i$ |
| $r$ | 0 for other features, 1 for this feature |
| $\alpha$ | baseline difference between this feature and the other features |
| $x_i$ | 0 for control, 1 for case |
| $\beta$ | condition-by-feature interaction |
| $x_i r$ | product that turns the interaction on only for case samples and this-feature rows |

This simplified model is the hand-expanded version of:

```r
~ sample + exon + condition:exon
```

The real DEXSeq model also has details for size factors, dispersions, multiple bins, constraints, and covariates. But this small version captures the key idea.

---

## 7. Why the sample term matters

The term:

$$
s_i
$$

is a sample-specific effect.

It lets each sample have its own overall gene expression level.

For example, sample C1 has gene total 100, while sample D1 has gene total 200.

DEXSeq does not put those totals into the model as fixed offsets. Instead, it estimates a sample effect that applies to both rows from the same sample:

| sample | row type | DEXSeq sample effect |
|---|---|---|
| C1 | this feature | $s_{C1}$ |
| C1 | other features | $s_{C1}$ |
| D1 | this feature | $s_{D1}$ |
| D1 | other features | $s_{D1}$ |

Because the same $s_i$ appears in both rows from the same sample, it captures sample-level gene abundance.

Then the interaction term asks whether the split between this feature and the other features changes by condition.

---

## 8. Plug in the row types

The model is:

$$
\log(\mu_{i r}) = s_i + \alpha r + \beta x_i r
$$

Now plug in the two possible row types.

### Other features

For the other-features row:

$$
r = 0
$$

So:

$$
\log(\mu_{i,other}) = s_i + \alpha \times 0 + \beta x_i \times 0
$$

Both terms involving $r$ disappear:

$$
\log(\mu_{i,other}) = s_i
$$

### This feature

For the this-feature row:

$$
r = 1
$$

So:

$$
\log(\mu_{i,this}) = s_i + \alpha \times 1 + \beta x_i \times 1
$$

which becomes:

$$
\log(\mu_{i,this}) = s_i + \alpha + \beta x_i
$$

So DEXSeq is fitting these two expected counts:

| row type | model |
|---|---|
| other features | $\log(\mu_{i,other}) = s_i$ |
| this feature | $\log(\mu_{i,this}) = s_i + \alpha + \beta x_i$ |

The same sample effect $s_i$ appears in both rows.

---

## 9. The sample effect cancels in the ratio

The key move is to subtract the other-features equation from the this-feature equation.

This-feature row:

$$
\log(\mu_{i,this}) = s_i + \alpha + \beta x_i
$$

Other-features row:

$$
\log(\mu_{i,other}) = s_i
$$

Subtract:

$$
\log(\mu_{i,this}) - \log(\mu_{i,other}) = s_i + \alpha + \beta x_i - s_i
$$

The $s_i$ terms cancel:

$$
\log(\mu_{i,this}) - \log(\mu_{i,other}) = \alpha + \beta x_i
$$

Use the log ratio rule:

$$
\log\left(\frac{\mu_{i,this}}{\mu_{i,other}}\right) = \alpha + \beta x_i
$$

This is the heart of DEXSeq.

The model says:

> log expected ratio of this feature to the rest of the gene = baseline ratio + condition interaction

So DEXSeq gets a usage-like comparison because the sample effect cancels out.

That is different from the offset model, where the gene total is subtracted as a fixed offset.

---

## 10. Work through F1

For F1, the observed average counts are:

| condition | F1 | other features |
|---|---:|---:|
| control | 80 | 20 |
| case | 80 | 120 |

The control ratio is:

$$
\frac{\text{F1}}{\text{other}} = \frac{80}{20} = 4
$$

The case ratio is:

$$
\frac{\text{F1}}{\text{other}} = \frac{80}{120} = \frac{2}{3}
$$

The DEXSeq interaction coefficient is the change in the log ratio:

$$
\beta_{F1,\ DEXSeq} = \log\left(\frac{2}{3}\right) - \log(4)
$$

Use the log ratio rule:

$$
\beta_{F1,\ DEXSeq} = \log\left(\frac{2/3}{4}\right)
$$

$$
\beta_{F1,\ DEXSeq} = \log\left(\frac{1}{6}\right)
$$

$$
\beta_{F1,\ DEXSeq} \approx -1.792
$$

Exponentiate:

$$
e^{\beta_{F1,\ DEXSeq}} = \frac{1}{6}
$$

Interpretation:

> The F1-to-other ratio in case is one-sixth of the F1-to-other ratio in control.

That is not the same number as the saseR-style offset coefficient:

$$
\beta_{F1,\ offset} = \log\left(\frac{0.4}{0.8}\right) = \log\left(\frac{1}{2}\right)
$$

Both say F1 usage decreased.

But they measure the effect on different scales:

| model | coefficient meaning | F1 value |
|---|---|---:|
| saseR-style offset | log usage ratio, case usage divided by control usage | $\log(1/2)$ |
| DEXSeq interaction | log odds ratio, case this-vs-other ratio divided by control this-vs-other ratio | $\log(1/6)$ |

---

## 11. Why DEXSeq's coefficient is an odds ratio

Feature usage is:

$$
p = \frac{\text{this feature}}{\text{gene total}}
$$

The rest of the gene has usage:

$$
1 - p
$$

DEXSeq compares this feature to the rest:

$$
\frac{\text{this feature}}{\text{other features}}
$$

In usage terms, that is:

$$
\frac{p}{1 - p}
$$

This is called an **odds**.

So DEXSeq's interaction coefficient is:

$$
\beta_{DEXSeq} = \log\left(\frac{p_{case}/(1-p_{case})}{p_{control}/(1-p_{control})}\right)
$$

That is a **log odds ratio**.

For F1:

$$
p_{control} = 0.8
$$

$$
p_{case} = 0.4
$$

Control odds:

$$
\frac{p_{control}}{1-p_{control}} = \frac{0.8}{0.2} = 4
$$

Case odds:

$$
\frac{p_{case}}{1-p_{case}} = \frac{0.4}{0.6} = \frac{2}{3}
$$

Odds ratio:

$$
\frac{2/3}{4} = \frac{1}{6}
$$

Log odds ratio:

$$
\log\left(\frac{1}{6}\right) \approx -1.792
$$

So DEXSeq is not directly estimating:

$$
\log\left(\frac{p_{case}}{p_{control}}\right)
$$

It is estimating a condition effect on:

$$
\log\left(\frac{p}{1-p}\right)
$$

for the feature being tested against the rest of the gene.

---

## 12. Work through F2

For F2, the observed average counts are:

| condition | F2 | other features |
|---|---:|---:|
| control | 20 | 80 |
| case | 120 | 80 |

The control ratio is:

$$
\frac{\text{F2}}{\text{other}} = \frac{20}{80} = \frac{1}{4}
$$

The case ratio is:

$$
\frac{\text{F2}}{\text{other}} = \frac{120}{80} = \frac{3}{2}
$$

The DEXSeq interaction coefficient is:

$$
\beta_{F2,\ DEXSeq} = \log\left(\frac{3/2}{1/4}\right)
$$

$$
\beta_{F2,\ DEXSeq} = \log(6)
$$

$$
\beta_{F2,\ DEXSeq} \approx 1.792
$$

Exponentiate:

$$
e^{\beta_{F2,\ DEXSeq}} = 6
$$

Interpretation:

> The F2-to-other ratio in case is 6 times the F2-to-other ratio in control.

Compare this with the offset coefficient:

$$
\beta_{F2,\ offset} = \log\left(\frac{0.6}{0.2}\right) = \log(3)
$$

Again, both models say F2 usage increased, but the coefficient is on a different scale.

---

## 13. The null hypothesis is the same idea

DEXSeq tests whether the interaction is zero.

The reduced model is:

```r
~ sample + exon
```

In the simplified algebra:

$$
\log(\mu_{i r}) = s_i + \alpha r
$$

This model says:

> The this-feature-to-other-features ratio is the same in control and case.

The full model is:

```r
~ sample + exon + condition:exon
```

In the simplified algebra:

$$
\log(\mu_{i r}) = s_i + \alpha r + \beta x_i r
$$

This model says:

> The this-feature-to-other-features ratio may differ between control and case.

The null hypothesis is:

$$
H_0: \beta = 0
$$

If:

$$
\beta = 0
$$

then:

$$
e^\beta = 1
$$

That means:

> The ratio of this feature to the rest of the gene is unchanged between conditions.

If the full model fits much better than the reduced model, DEXSeq reports evidence for differential usage.

---

## 14. Why an interaction is needed

To see why the interaction is needed, compare three models.

### Model A: sample only

```r
~ sample
```

This allows every sample to have its own expression level.

But it does not distinguish this feature from the other features.

### Model B: sample plus exon

```r
~ sample + exon
```

This allows this feature to have a different baseline count from the other features.

But it assumes that difference is the same in control and case.

In ratio language:

$$
\log\left(\frac{\mu_{i,this}}{\mu_{i,other}}\right) = \alpha
$$

There is no $x_i$ in that equation.

So condition cannot change the ratio.

### Model C: sample plus exon plus condition-by-exon interaction

```r
~ sample + exon + condition:exon
```

This allows the this-feature-versus-other-features difference to change by condition.

In ratio language:

$$
\log\left(\frac{\mu_{i,this}}{\mu_{i,other}}\right) = \alpha + \beta x_i
$$

Now condition can change the ratio.

That is why the DTU term is the interaction term.

---

## 15. DEXSeq paper notation

The DEXSeq paper writes the model with more indices because it models genes, samples, and counting bins.

The paper uses:

$$
K_{i j l} = \text{count for gene } i \text{, sample } j \text{, counting bin } l
$$

and a negative-binomial count model:

$$
K_{i j l} \sim NB(\text{mean} = s_j \mu_{i j l}, \text{dispersion} = \alpha_{i l})
$$

Here:

| term | meaning |
|---|---|
| $i$ | gene |
| $j$ | sample |
| $l$ | counting bin, such as an exon bin |
| $s_j$ | sequencing-depth size factor for sample $j$ |
| $\mu_{i j l}$ | normalized expected abundance for the bin |
| $\alpha_{i l}$ | dispersion for this gene-bin |

The DEXSeq paper's usage model can be written as:

$$
\log(\mu_{i j l}) = \beta_i^G + \beta_{i l}^E + \beta_{i j}^S + \beta_{i r_j l}^{EC}
$$

The terms mean:

| term | rough meaning |
|---|---|
| $\beta_i^G$ | baseline expression strength of gene $i$ |
| $\beta_{i l}^E$ | baseline effect for counting bin $l$ |
| $\beta_{i j}^S$ | sample-specific effect |
| $\beta_{i r_j l}^{EC}$ | exon/bin-by-condition interaction |

The most important term for DTU is:

$$
\beta_{i r_j l}^{EC}
$$

That is the interaction term.

If this interaction is zero, then condition does not change the usage of bin $l$.

If this interaction is nonzero, then condition changes the usage of bin $l$ relative to the gene.

---

## 16. DEXSeq versus the saseR-style offset model

The two models can be lined up like this.

### saseR-style offset model

For one feature:

$$
\log(\mu_{f i}) = \log(G_i) + \alpha_f + \beta_f x_i
$$

Subtract the offset:

$$
\log\left(\frac{\mu_{f i}}{G_i}\right) = \alpha_f + \beta_f x_i
$$

So the condition coefficient is about:

$$
\frac{\text{feature}}{\text{gene total}}
$$

That is usage.

### DEXSeq interaction model

For one feature-versus-other comparison:

$$
\log(\mu_{i,this}) = s_i + \alpha + \beta x_i
$$

$$
\log(\mu_{i,other}) = s_i
$$

Subtract:

$$
\log\left(\frac{\mu_{i,this}}{\mu_{i,other}}\right) = \alpha + \beta x_i
$$

So the condition coefficient is about:

$$
\frac{\text{feature}}{\text{other features in the gene}}
$$

That is an odds-like usage ratio.

The core difference is:

| model | normalization idea | coefficient tests |
|---|---|---|
| saseR-style offset | divide by gene total using a fixed offset | feature usage relative to gene total |
| DEXSeq interaction | compare feature count to other-gene count using sample blocking | feature-to-other ratio |

---

## 17. Why DEXSeq can be slower

The DEXSeq model includes a sample-specific term:

$$
s_i
$$

In R formula form, that is the:

```r
sample
```

part of:

```r
~ sample + exon + condition:exon
```

That means the model matrix grows as the number of samples grows.

If there are 6 samples, there is roughly one sample-level coefficient per sample, subject to the usual model constraints that prevent duplicate intercepts.

If there are 1,000 samples, there are roughly 1,000 sample-level coefficients.

Those sample terms are useful because they absorb sample-level gene expression differences. But they also make the model larger.

The saseR-style offset model avoids estimating a large set of sample blocking coefficients in each gene model. It puts the gene total into the model as a known offset:

$$
\log(G_i)
$$

So instead of estimating a sample effect and letting it cancel in a ratio, the offset model directly subtracts the gene total on the log scale.

This is the scalability idea emphasized in the saseR paper.

---

## 18. A simple analogy

Imagine a gene as a pie.

The gene total is the size of the whole pie.

The feature count is the size of one slice.

### Offset model

The offset model asks:

> How big is this slice compared with the whole pie?

Mathematically:

$$
\frac{\text{slice}}{\text{whole pie}}
$$

### DEXSeq model

DEXSeq asks:

> How big is this slice compared with the rest of the pie?

Mathematically:

$$
\frac{\text{slice}}{\text{rest of pie}}
$$

Both can detect that the slice changed.

But the reported effect size is different because the denominator is different.

---

## 19. Important interpretation warning

Do not interpret a DEXSeq interaction coefficient as:

$$
\log\left(\frac{p_{case}}{p_{control}}\right)
$$

That is the offset-model interpretation.

For the simplified DEXSeq feature-versus-other setup, the coefficient is better thought of as:

$$
\log\left(\frac{p_{case}/(1-p_{case})}{p_{control}/(1-p_{control})}\right)
$$

That is:

> log odds ratio for this feature versus the rest of the gene

The p-value still answers the DTU question:

> Is there evidence that usage changed?

But the coefficient is not on the same scale as the saseR-style offset coefficient.

---

## 20. Summary

DEXSeq and the saseR-style offset model are both GLM approaches for usage.

The saseR-style offset model says:

$$
\log(\mu_{f i}) = \log(G_i) + \alpha_f + \beta_f x_i
$$

This directly models:

$$
\log\left(\frac{\mu_{f i}}{G_i}\right)
$$

So $\beta_f$ is a log usage-ratio effect.

DEXSeq says, in simplified form:

$$
\log(\mu_{i r}) = s_i + \alpha r + \beta x_i r
$$

This implies:

$$
\log\left(\frac{\mu_{i,this}}{\mu_{i,other}}\right) = \alpha + \beta x_i
$$

So $\beta$ is a condition-by-feature interaction, interpretable as a log odds-ratio effect for this feature versus the rest of the gene.

The big conceptual difference is:

| question | saseR-style offset | DEXSeq |
|---|---|---|
| How is gene expression handled? | fixed gene-total offset | estimated sample blocking term |
| What is compared? | feature versus gene total | feature versus other features |
| Key term | condition coefficient $\beta_f$ | condition-by-exon interaction $\beta$ |
| Effect scale | log usage ratio | log odds ratio |
| Reduced model | no condition effect on feature usage | no condition-by-exon interaction |
| Full model | condition effect allowed | condition-by-exon interaction allowed |

The key DEXSeq idea is:

> The interaction term asks whether condition changes the balance between this feature and the rest of the gene.
