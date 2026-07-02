# GLMs for Differential Transcript Usage with a Gene-Expression Offset

This note walks through a toy example of how generalized linear models (GLMs) can be used for **differential transcript usage (DTU)**, including the role of a **gene-expression offset**. The same logic applies to transcript isoforms, exon bins, splice junctions, splice-graph features, or DSF-style splicing features.

The central idea is:

> **Ordinary transcript differential expression:** “Did this transcript’s raw count change?”  
> **DTU with a gene-expression offset:** “Did this transcript’s count change more or less than expected, given the parent gene’s total expression?”

In tools such as **saseR**, this kind of offset logic is used to adapt standard negative-binomial RNA-seq GLMs, such as those from edgeR or DESeq2, to test **feature usage relative to the parent gene** rather than raw feature abundance.

---

## 1. The toy data

Suppose we have one gene with two transcript isoforms:

- Transcript **T1**
- Transcript **T2**

We compare two groups:

- **Control**
- **Case**

Here are six samples:

| sample | condition | gene total counts | T1 counts | T2 counts |
|---|---|---:|---:|---:|
| C1 | control | 100 | 80 | 20 |
| C2 | control | 120 | 96 | 24 |
| C3 | control | 80 | 64 | 16 |
| D1 | case | 200 | 80 | 120 |
| D2 | case | 180 | 72 | 108 |
| D3 | case | 220 | 88 | 132 |

The gene total is:

$$
G_i = Y_{T1,i} + Y_{T2,i}
$$

For sample C1:

$$
G_{C1} = 80 + 20 = 100
$$

For sample D1:

$$
G_{D1} = 80 + 120 = 200
$$

So the gene itself is more highly expressed in cases.

Average counts:

| condition | gene total | T1 counts | T2 counts |
|---|---:|---:|---:|
| control | 100 | 80 | 20 |
| case | 200 | 80 | 120 |

A normal differential expression test on transcript counts would notice:

$$
T1: 80 \rightarrow 80
$$

so T1 looks unchanged.

For T2:

$$
T2: 20 \rightarrow 120
$$

so T2 looks strongly increased.

But that is not the DTU question.

---

## 2. Convert counts to transcript usage

Transcript usage means:

$$
\text{transcript usage} = \frac{\text{transcript count}}{\text{total gene count}}
$$

For transcript T1 in control:

$$
\text{usage}_{T1,\ control} = \frac{80}{100} = 0.8
$$

So T1 is 80% of the gene's expression in control.

For transcript T1 in case:

$$
\text{usage}_{T1,\ case} = \frac{80}{200} = 0.4
$$

So T1 is only 40% of the gene's expression in case.

For T2:

$$
\text{usage}_{T2,\ control} = \frac{20}{100} = 0.2
$$

$$
\text{usage}_{T2,\ case} = \frac{120}{200} = 0.6
$$

So the gene shifts from mostly T1 to mostly T2.

| transcript | control usage | case usage |
|---|---:|---:|
| T1 | 80% | 40% |
| T2 | 20% | 60% |

This is the key point:

| transcript | raw-count interpretation | usage interpretation |
|---|---|---|
| T1 | unchanged | usage decreased |
| T2 | increased | usage increased |

A DTU model wants the second interpretation.

---

## 3. Define every term

We will focus on one transcript at a time.

Let:

$$
i = \text{sample index}
$$

So:

$$
i = 1,2,3,4,5,6
$$

Let:

$$
t = \text{transcript index}
$$

Here:

$$
t \in \{T1, T2\}
$$

Let:

$$
Y_{t i} = \text{observed read count for transcript } t \text{ in sample } i
$$

For example:

$$
Y_{T1, C1} = 80
$$

$$
Y_{T2, C1} = 20
$$

Let:

$$
G_i = \text{total expression count for the parent gene in sample } i
$$

For C1:

$$
G_{C1} = 100
$$

For D1:

$$
G_{D1} = 200
$$

Let:

$$
x_i = \text{condition indicator for sample } i
$$

We encode condition as:

$$
x_i =
\begin{cases}
0 & \text{if sample } i \text{ is control} \\
1 & \text{if sample } i \text{ is case}
\end{cases}
$$

So:

| sample | condition | $x_i$ |
|---|---|---:|
| C1 | control | 0 |
| C2 | control | 0 |
| C3 | control | 0 |
| D1 | case | 1 |
| D2 | case | 1 |
| D3 | case | 1 |

Let:

$$
\mu_{t i} = \text{expected count for transcript } t \text{ in sample } i
$$

This is not the observed count. It is what the model predicts on average.

Let:

$$
p_{t i} = \text{expected usage of transcript } t \text{ in sample } i
$$

So, approximately:

$$
p_{t i} \approx \frac{Y_{t i}}{G_i}
$$

The core relationship is:

$$
\mu_{t i} = G_i \times p_{t i}
$$

In words:

> expected transcript count = total gene expression × transcript usage

That is the whole idea behind the offset model.

---

## 4. Start without a GLM

For T1 in control:

$$
G_i = 100
$$

$$
p_{T1,i} = 0.8
$$

Therefore:

$$
\mu_{T1,i} = G_i \times p_{T1,i}
$$

$$
\mu_{T1,i} = 100 \times 0.8
$$

$$
\mu_{T1,i} = 80
$$

For T1 in case:

$$
G_i = 200
$$

$$
p_{T1,i} = 0.4
$$

$$
\mu_{T1,i} = 200 \times 0.4
$$

$$
\mu_{T1,i} = 80
$$

So T1 raw counts can be unchanged even though T1 usage decreased.

That is why the offset is needed.

---

## 5. Why use logs?

RNA-seq GLMs usually model counts using a **log link**.

Instead of modeling:

$$
\mu_{t i}
$$

directly, we model:

$$
\log(\mu_{t i})
$$

Why? Because counts must be positive, and exponentiating a number always gives a positive value.

For example:

$$
e^0 = 1
$$

$$
e^1 \approx 2.718
$$

$$
e^{-1} \approx 0.368
$$

So if the model predicts something on the log scale, we can convert it back to the count scale using:

$$
\mu_{t i} = e^{\text{linear predictor}}
$$

---

## 6. Deriving the offset model

We start from:

$$
\mu_{t i} = G_i \times p_{t i}
$$

Take the log of both sides:

$$
\log(\mu_{t i}) = \log(G_i \times p_{t i})
$$

A log rule says:

$$
\log(a \times b) = \log(a) + \log(b)
$$

So:

$$
\log(\mu_{t i}) = \log(G_i) + \log(p_{t i})
$$

So far this is just algebra. The expected count splits cleanly into two pieces:

- $\log(G_i)$, the part coming from how much the **whole gene** was expressed in sample $i$, and
- $\log(p_{t i})$, the part coming from this **transcript's usage** in sample $i$.

We have not yet said anything about *what controls* usage. That is the next step, and it is where the modeling decision happens.

### Step 1: We still need a formula for $\log(p_{t i})$

Right now $\log(p_{t i})$ is just a placeholder. It is a different number for every sample, and nothing ties those numbers together. To actually fit a model and run a test, we need to express $p_{t i}$ in terms of something we know about each sample.

The only thing that distinguishes our samples is their **condition** (control vs. case), encoded by $x_i$:

$$
x_i =
\begin{cases}
0 & \text{control} \\
1 & \text{case}
\end{cases}
$$

So we want usage to be a function of $x_i$.

### Step 2: Why model usage on the log scale

We choose to write the model for $\log(p_{t i})$ rather than for $p_{t i}$ directly. The reason is the same as for the count mean: a log scale turns **multiplicative** effects into **additive** ones.

A "usage doubled" or "usage halved" statement is multiplicative on the raw scale:

$$
p_{t i} = p_{\text{baseline}} \times (\text{fold change})
$$

Taking logs makes it additive:

$$
\log(p_{t i}) = \log(p_{\text{baseline}}) + \log(\text{fold change})
$$

Additive effects are exactly what a linear model can represent.

### Step 3: The simplest log-linear model

The simplest function of $x_i$ is a straight line:

$$
\log(p_{t i}) = \alpha_t + \beta_t x_i
$$

Here the two unknowns mean:

- $\alpha_t$ = the **intercept**: the value of $\log(p_{t i})$ when $x_i = 0$, i.e. baseline log usage of transcript $t$ in controls.
- $\beta_t$ = the **slope**: how much $\log(p_{t i})$ changes when we go from control ($x_i = 0$) to case ($x_i = 1$).

We can see this by plugging in each condition:

$$
\text{control: } \log(p_{t i}) = \alpha_t + \beta_t \times 0 = \alpha_t
$$

$$
\text{case: } \log(p_{t i}) = \alpha_t + \beta_t \times 1 = \alpha_t + \beta_t
$$

So the difference between the two conditions is exactly $\beta_t$:

$$
\beta_t = \log(p_{t,\ case}) - \log(p_{t,\ control}) = \log\!\left(\frac{p_{t,\ case}}{p_{t,\ control}}\right)
$$

In words, $\beta_t$ is the **log fold change in usage** between case and control. That single number is the thing the DTU test ultimately cares about.

### Step 4: Substitute the usage model back into the count model

Now we have two equations. The first came from pure algebra:

$$
\log(\mu_{t i}) = \log(G_i) + \log(p_{t i})
$$

The second is our chosen model for usage:

$$
\log(p_{t i}) = \alpha_t + \beta_t x_i
$$

We substitute the second into the first by replacing $\log(p_{t i})$ with $\alpha_t + \beta_t x_i$:

$$
\log(\mu_{t i}) = \log(G_i) + \underbrace{\big(\alpha_t + \beta_t x_i\big)}_{\log(p_{t i})}
$$

which gives:

$$
\log(\mu_{t i}) = \log(G_i) + \alpha_t + \beta_t x_i
$$

This is the DTU offset GLM. The key thing to notice is that the gene-expression term $\log(G_i)$ came in for free from the algebra in Step 0, while $\alpha_t$ and $\beta_t$ are the only quantities the model actually estimates.

---

## 7. What each term means

The model is:

$$
\log(\mu_{t i}) = \log(G_i) + \alpha_t + \beta_t x_i
$$

Term by term:

| term | meaning |
|---|---|
| $Y_{t i}$ | observed count for transcript $t$ in sample $i$ |
| $\mu_{t i}$ | expected count for transcript $t$ in sample $i$ |
| $G_i$ | total parent-gene expression in sample $i$ |
| $\log(G_i)$ | the offset |
| $\alpha_t$ | baseline log usage of transcript $t$ in controls |
| $\beta_t$ | change in log usage of transcript $t$ in cases versus controls |
| $x_i$ | condition indicator: 0 for control, 1 for case |
| $\log$ | natural logarithm |
| $e^z$ | exponential function, the inverse of $\log(z)$ |

The offset is:

$$
\log(G_i)
$$

It is called an **offset** because its coefficient is fixed at 1.

That means the model is really:

$$
\log(\mu_{t i}) = 1 \times \log(G_i) + \alpha_t + \beta_t x_i
$$

The model is not allowed to estimate the coefficient on $\log(G_i)$. It is forced to be 1.

That encodes the idea:

> if gene expression doubles, expected transcript counts should double, unless usage also changes.

---

## 8. Work through T1 fully

For T1, the usage values are:

| condition | T1 usage |
|---|---:|
| control | 0.8 |
| case | 0.4 |

The model says:

$$
\log(p_{T1,i}) = \alpha_{T1} + \beta_{T1} x_i
$$

### Controls

For controls:

$$
x_i = 0
$$

So:

$$
\log(p_{T1,control}) = \alpha_{T1} + \beta_{T1} \times 0
$$

$$
\log(p_{T1,control}) = \alpha_{T1}
$$

But we know:

$$
p_{T1,control} = 0.8
$$

So:

$$
\alpha_{T1} = \log(0.8)
$$

$$
\alpha_{T1} \approx -0.223
$$

So the baseline log usage for T1 is:

$$
\alpha_{T1} \approx -0.223
$$

### Cases

For cases:

$$
x_i = 1
$$

So:

$$
\log(p_{T1,case}) = \alpha_{T1} + \beta_{T1} \times 1
$$

$$
\log(p_{T1,case}) = \alpha_{T1} + \beta_{T1}
$$

We know:

$$
p_{T1,case} = 0.4
$$

So:

$$
\log(0.4) = \alpha_{T1} + \beta_{T1}
$$

Substitute:

$$
\alpha_{T1} = \log(0.8)
$$

$$
\log(0.4) = \log(0.8) + \beta_{T1}
$$

Now solve for $\beta_{T1}$:

$$
\beta_{T1} = \log(0.4) - \log(0.8)
$$

Another log rule says:

$$
\log(a) - \log(b) = \log\left(\frac{a}{b}\right)
$$

So:

$$
\beta_{T1} = \log\left(\frac{0.4}{0.8}\right)
$$

$$
\beta_{T1} = \log(0.5)
$$

$$
\beta_{T1} \approx -0.693
$$

So the model estimates:

$$
\beta_{T1} \approx -0.693
$$

On the natural scale:

$$
e^{\beta_{T1}} = e^{-0.693} \approx 0.5
$$

Why does exponentiating $\beta_{T1}$ give a usage ratio? Recall from Step 3 that

$$
\beta_{T1} = \log\!\left(\frac{p_{T1,\ case}}{p_{T1,\ control}}\right)
$$

so $\beta_{T1}$ is already the **log** of the ratio of the two usage values. Exponentiating simply undoes that log and returns the ratio itself:

$$
e^{\beta_{T1}} = e^{\log\left(\frac{p_{T1,\ case}}{p_{T1,\ control}}\right)} = \frac{p_{T1,\ case}}{p_{T1,\ control}}
$$

We can check this against the usage symbols defined earlier. From Section 2:

$$
p_{T1,\ control} = 0.8, \qquad p_{T1,\ case} = 0.4
$$

Substituting:

$$
e^{\beta_{T1}} = \frac{p_{T1,\ case}}{p_{T1,\ control}} = \frac{0.4}{0.8} = 0.5
$$

which matches the value above.

That means:

> T1 usage in case is 0.5 times its usage in control.

In other words, $p_{T1,\ case}$ is half of $p_{T1,\ control}$ — T1 usage is cut in half.

---

## 9. Use the fitted T1 model to predict counts

The model is:

$$
\log(\mu_{T1,i}) = \log(G_i) + \alpha_{T1} + \beta_{T1} x_i
$$

We found:

$$
\alpha_{T1} = \log(0.8)
$$

$$
\beta_{T1} = \log(0.5)
$$

### Predict T1 count for a control sample with $G_i = 100$

Because it is control:

$$
x_i = 0
$$

$$
\log(\mu_{T1,i}) = \log(100) + \log(0.8) + \log(0.5) \times 0
$$

The last term is zero:

$$
\log(0.5) \times 0 = 0
$$

So:

$$
\log(\mu_{T1,i}) = \log(100) + \log(0.8)
$$

Using the log product rule in reverse:

$$
\log(100) + \log(0.8) = \log(100 \times 0.8)
$$

$$
\log(\mu_{T1,i}) = \log(80)
$$

Therefore:

$$
\mu_{T1,i} = 80
$$

### Predict T1 count for a case sample with $G_i = 200$

Because it is case:

$$
x_i = 1
$$

$$
\log(\mu_{T1,i}) = \log(200) + \log(0.8) + \log(0.5) \times 1
$$

$$
\log(\mu_{T1,i}) = \log(200) + \log(0.8) + \log(0.5)
$$

Combine the logs:

$$
\log(\mu_{T1,i}) = \log(200 \times 0.8 \times 0.5)
$$

$$
\log(\mu_{T1,i}) = \log(80)
$$

Therefore:

$$
\mu_{T1,i} = 80
$$

So the GLM says:

$$
\text{case T1 count} = 200 \times 0.8 \times 0.5 = 80
$$

Interpretation:

- The gene doubled: $100 \rightarrow 200$
- T1 usage halved: $0.8 \rightarrow 0.4$
- These cancel out
- So the raw T1 count stays at 80

That is exactly why raw transcript DE and DTU can disagree.

---

## 10. Work through T2 fully

For T2:

| condition | T2 usage |
|---|---:|
| control | 0.2 |
| case | 0.6 |

The same model is:

$$
\log(p_{T2,i}) = \alpha_{T2} + \beta_{T2} x_i
$$

### Controls

For controls:

$$
x_i = 0
$$

$$
\log(p_{T2,control}) = \alpha_{T2}
$$

Since:

$$
p_{T2,control} = 0.2
$$

we get:

$$
\alpha_{T2} = \log(0.2)
$$

$$
\alpha_{T2} \approx -1.609
$$

### Cases

For cases:

$$
x_i = 1
$$

$$
\log(p_{T2,case}) = \alpha_{T2} + \beta_{T2}
$$

Since:

$$
p_{T2,case} = 0.6
$$

we get:

$$
\log(0.6) = \log(0.2) + \beta_{T2}
$$

Solve:

$$
\beta_{T2} = \log(0.6) - \log(0.2)
$$

$$
\beta_{T2} = \log\left(\frac{0.6}{0.2}\right)
$$

$$
\beta_{T2} = \log(3)
$$

$$
\beta_{T2} \approx 1.099
$$

So:

$$
e^{\beta_{T2}} = e^{1.099} \approx 3
$$

Interpretation:

> T2 usage is 3 times higher in case than in control.

---

## 11. Why T2 raw counts increased 6-fold

T2 raw count changed from:

$$
20 \rightarrow 120
$$

That is:

$$
\frac{120}{20} = 6
$$

So raw T2 expression increased 6-fold.

But the GLM decomposes that into two parts:

$$
\text{raw transcript change} = \text{gene expression change} \times \text{usage change}
$$

Gene expression changed from 100 to 200:

$$
\frac{200}{100} = 2
$$

T2 usage changed from 0.2 to 0.6:

$$
\frac{0.6}{0.2} = 3
$$

So:

$$
6 = 2 \times 3
$$

The offset accounts for the $2\times$ gene expression increase. The coefficient $\beta_{T2}$ captures the remaining $3\times$ usage increase.

---

## 12. What would happen without the offset?

Without the offset, the model would be:

$$
\log(\mu_{t i}) = \alpha_t + \beta_t x_i
$$

This asks:

> Did the raw transcript count change between conditions?

For T2:

$$
\beta_{T2,\ raw} = \log\left(\frac{120}{20}\right)
$$

$$
\beta_{T2,\ raw} = \log(6)
$$

$$
\beta_{T2,\ raw} \approx 1.792
$$

With the offset, the model asks:

> Did transcript usage change after accounting for total gene expression?

For T2:

$$
\beta_{T2,\ usage} = \log\left(\frac{0.6}{0.2}\right)
$$

$$
\beta_{T2,\ usage} = \log(3)
$$

$$
\beta_{T2,\ usage} \approx 1.099
$$

So:

| model | question | T2 coefficient |
|---|---|---:|
| no offset | raw transcript expression changed? | $\log(6)$ |
| gene offset | transcript usage changed? | $\log(3)$ |

The offset removes the part explained by whole-gene expression.

---

## 13. The GLM as a count model

So far we talked about expected counts. A GLM also needs a probability model for the observed counts.

The common RNA-seq model is the negative binomial:

$$
Y_{t i} \sim NB(\mu_{t i}, \phi_t)
$$

This means:

> observed transcript counts $Y_{t i}$ are treated as noisy observations around expected counts $\mu_{t i}$.

The terms are:

| term | meaning |
|---|---|
| $Y_{t i}$ | observed count |
| $NB$ | negative-binomial distribution |
| $\mu_{t i}$ | expected count |
| $\phi_t$ | dispersion for transcript $t$ |
| $\sim$ | “is distributed as” |

The negative binomial is used because RNA-seq counts are usually more variable than a simple Poisson model expects.

The variance is often written as:

$$
Var(Y_{t i}) = \mu_{t i} + \phi_t \mu_{t i}^2
$$

Term by term:

| term | meaning |
|---|---|
| $Var(Y_{t i})$ | variance of the observed count |
| $\mu_{t i}$ | expected count |
| $\phi_t$ | dispersion |
| $\phi_t \mu_{t i}^2$ | extra biological/technical variability beyond Poisson noise |

If:

$$
\phi_t = 0
$$

then:

$$
Var(Y_{t i}) = \mu_{t i}
$$

That is the Poisson case.

If:

$$
\phi_t > 0
$$

then:

$$
Var(Y_{t i}) > \mu_{t i}
$$

That is overdispersion, which is common in RNA-seq.

---

## 14. The full DTU GLM

Putting the pieces together:

$$
Y_{t i} \sim NB(\mu_{t i}, \phi_t)
$$

$$
\log(\mu_{t i}) = \log(G_i) + \alpha_t + \beta_t x_i
$$

This is the full model.

It has two layers.

### Count layer

$$
Y_{t i} \sim NB(\mu_{t i}, \phi_t)
$$

This says:

> observed counts are noisy negative-binomial measurements.

### Mean layer

$$
\log(\mu_{t i}) = \log(G_i) + \alpha_t + \beta_t x_i
$$

This says:

> expected transcript counts depend on total gene expression plus a condition effect on usage.

---

## 15. Where the hypothesis test comes in

For each transcript, the DTU test asks:

$$
H_0: \beta_t = 0
$$

This is the null hypothesis.

Recall from Step 3 what $\beta_t$ means in terms of usage:

$$
\beta_t = \log\!\left(\frac{p_{t,\ case}}{p_{t,\ control}}\right)
$$

So $\beta_t$ is the **log fold change in usage** between case and control. Setting it to zero,

$$
\beta_t = 0 \quad\Longleftrightarrow\quad \log\!\left(\frac{p_{t,\ case}}{p_{t,\ control}}\right) = 0 \quad\Longleftrightarrow\quad \frac{p_{t,\ case}}{p_{t,\ control}} = 1 \quad\Longleftrightarrow\quad p_{t,\ case} = p_{t,\ control}
$$

means the two usage values are identical.

It means:

> transcript usage does not differ between case and control.

Because:

$$
e^{\beta_t} = e^0 = 1
$$

So if:

$$
\beta_t = 0
$$

then:

$$
e^{\beta_t} = 1
$$

and case usage is 1 times control usage, meaning no usage change.

The alternative hypothesis is:

$$
H_A: \beta_t \ne 0
$$

This means:

> transcript usage differs between case and control.

For T1:

$$
\beta_{T1} = \log(0.5)
$$

So:

$$
\beta_{T1} < 0
$$

T1 usage decreases.

For T2:

$$
\beta_{T2} = \log(3)
$$

So:

$$
\beta_{T2} > 0
$$

T2 usage increases.

---

## 16. A likelihood-ratio version of the test

A GLM can compare two models.

### Reduced model

The reduced model has no condition effect:

$$
\log(\mu_{t i}) = \log(G_i) + \alpha_t
$$

This model says:

> transcript usage is the same in control and case.

### Full model

The full model includes condition:

$$
\log(\mu_{t i}) = \log(G_i) + \alpha_t + \beta_t x_i
$$

This model says:

> transcript usage may differ between control and case.

Then the software compares how well the two models fit.

The likelihood-ratio statistic is:

$$
LRT = 2 \times \left[\log L(\text{full model}) - \log L(\text{reduced model})\right]
$$

Term by term:

| term | meaning |
|---|---|
| $LRT$ | likelihood-ratio test statistic |
| $L$ | likelihood, meaning how plausible the observed data are under a model |
| $\log L$ | log likelihood |
| full model | model with condition effect |
| reduced model | model without condition effect |

If the full model fits much better, then:

$$
\log L(\text{full model}) \gg \log L(\text{reduced model})
$$

and the LRT statistic becomes large.

A large LRT statistic gives a small p-value.

A small p-value means:

> allowing transcript usage to differ by condition explains the data much better.

---

## 17. Why the offset makes this a usage test

Return to the key equation:

$$
\log(\mu_{t i}) = \log(G_i) + \alpha_t + \beta_t x_i
$$

Subtract the offset from both sides:

$$
\log(\mu_{t i}) - \log(G_i) = \alpha_t + \beta_t x_i
$$

Use the log ratio rule:

$$
\log(a) - \log(b) = \log\left(\frac{a}{b}\right)
$$

So:

$$
\log\left(\frac{\mu_{t i}}{G_i}\right) = \alpha_t + \beta_t x_i
$$

But:

$$
\frac{\mu_{t i}}{G_i}
$$

is the expected transcript fraction relative to the gene.

So the model is really:

$$
\log(\text{expected transcript usage}) = \alpha_t + \beta_t x_i
$$

That is the most important derivation.

The offset version:

$$
\log(\mu_{t i}) = \log(G_i) + \alpha_t + \beta_t x_i
$$

is mathematically equivalent to:

$$
\log\left(\frac{\mu_{t i}}{G_i}\right) = \alpha_t + \beta_t x_i
$$

So the offset turns a count model into a usage model.

---

## 18. Why not just divide counts by gene expression and run a t-test?

You could compute:

$$
\frac{Y_{t i}}{G_i}
$$

and compare proportions directly.

But the GLM has advantages:

1. It models count uncertainty.
2. It handles low counts better.
3. It can include dispersion.
4. It can include covariates.
5. It works naturally with edgeR/DESeq2-style RNA-seq frameworks.

For example, a richer model could be:

$$
\log(\mu_{t i}) =
\log(G_i)
+
\alpha_t
+
\beta_t x_i
+
\gamma_t b_i
$$

where:

| term | meaning |
|---|---|
| $b_i$ | batch indicator or batch covariate |
| $\gamma_t$ | transcript-specific batch effect |

Then:

$$
\beta_t
$$

still tests condition-associated usage change, after accounting for both gene expression and batch.

---

## 19. How this generalizes to exon bins, junctions, or splice-graph features

The same idea works for any feature inside a gene.

The feature could be:

| feature type | count $Y_{f i}$ means |
|---|---|
| transcript isoform | reads assigned to transcript |
| exon bin | reads overlapping exon bin |
| splice junction | reads supporting junction |
| splice-graph segment | reads supporting graph feature |
| DSF-style feature | reads supporting differential splicing feature |

Then use:

$$
Y_{f i} = \text{count for feature } f \text{ in sample } i
$$

$$
G_i = \text{total expression of the parent gene in sample } i
$$

The GLM becomes:

$$
Y_{f i} \sim NB(\mu_{f i}, \phi_f)
$$

$$
\log(\mu_{f i}) = \log(G_i) + \alpha_f + \beta_f x_i
$$

Same interpretation:

$$
\beta_f
$$

is the condition effect on **feature usage relative to the gene**.

---

## 20. Intuitive summary

The offset says:

> First, tell the model how much expression the gene had. Then ask whether this feature got more or fewer reads than expected from that gene expression.

For T1:

| quantity | control | case |
|---|---:|---:|
| gene expression | 100 | 200 |
| T1 raw count | 80 | 80 |
| T1 usage | 80% | 40% |

Without the offset, T1 looks unchanged.

With the offset, T1 is strongly depleted relative to the gene.

For T2:

| quantity | control | case |
|---|---:|---:|
| gene expression | 100 | 200 |
| T2 raw count | 20 | 120 |
| T2 usage | 20% | 60% |

Without the offset, T2 has a 6-fold raw increase.

With the offset, T2 has a 3-fold usage increase.

So the offset removes the gene-level expression change and leaves the isoform-usage change.

The one-line version is:

$$
\log(\mu_{t i}) = \log(G_i) + \alpha_t + \beta_t x_i
$$

which is equivalent to:

$$
\log\left(\frac{\mu_{t i}}{G_i}\right) = \alpha_t + \beta_t x_i
$$

That second form makes the meaning clear:

> the GLM is modeling transcript usage, not raw transcript abundance.
