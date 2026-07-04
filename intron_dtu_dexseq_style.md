# Intron Usage Testing: From Offset GLMs to a DEXSeq-Style Paired-Count Framework

This note examines how the offset-based intron splicing models used in
**Diff-Splice-Finder (DSF)** can be reframed as a **DEXSeq-style paired-count
negative-binomial (NB) model**. The goal is to make the statistical assumptions
explicit, clarify what each denominator choice implies for the "other" count, and
identify where the two formulations agree and where they differ.

---

## 1. The shared conceptual structure

Both the DTU offset model (as in saseR / edgeR) and the DEXSeq paired-count model
decompose the observed count for a focal feature into two parts:

$$
\text{total available evidence} = \text{focal feature count} + \text{other count}
$$

The **offset approach** treats the total as a fixed, observed denominator $D_{i,s}$
(one number per intron per sample) and conditions the model on it:

$$
\log(\mu_{i,s}) = \beta_0^{(i)} + \beta_1^{(i)} \cdot \text{Group}_s + \log(D_{i,s})
$$

The **DEXSeq approach** treats both the focal count and the remainder count as
separate NB random variables and models them jointly:

| row | count | NB mean |
|---|---|---|
| focal | $Y_{i,s}$ | $\mu_{i,s}$ |
| other | $Y_{\text{other},i,s} = D_{i,s} - Y_{i,s}$ | $\mu_{\text{other},i,s}$ |

The choice of $D_{i,s}$ — what the "total available evidence" is — is the critical
modeling decision, and DSF has developed several answers to that question. Each one
implies a different biological meaning for $Y_{\text{other}}$.

---

## 2. DSF denominator modes and what they imply for $Y_{\text{other}}$

### 2a. `cluster_max` denominator

**Definition**

Introns that share a donor (5' splice site) form a *donor cluster*; introns that
share an acceptor (3' splice site) form an *acceptor cluster*.

$$
D_{i,s}^{\text{cluster}} = \max\!\left(\sum_{j \in \text{donor cluster}(i)} Y_{j,s},\; \sum_{j \in \text{acceptor cluster}(i)} Y_{j,s}\right)
$$

The maximum is taken because novel splice sites can be singletons in one cluster
dimension but participate in a well-populated cluster in the other; using the max
prevents a singleton denominator of 1 or 2 reads from inflating PSI to near 1 in
both conditions and masking a real change (see the *singleton artifact* discussion
below).

**The implied $Y_{\text{other}}$**

$$
Y_{\text{other},i,s}^{\text{cluster}} = D_{i,s}^{\text{cluster}} - Y_{i,s}
$$

This "other" count represents *all other introns competing at the same splice site
(donor or acceptor, whichever has higher total expression)*. In the DEXSeq framing,
both $Y_{i,s}$ and $Y_{\text{other},i,s}$ are observed junction counts, so both are
naturally NB-distributed with the same splice-site expression as their shared scale
factor.

**Biological interpretation of $\beta_1$**

$$
e^{\beta_1} \approx \frac{p_{i,\text{case}}}{p_{i,\text{control}}}, \quad
p_{i,s} = \frac{Y_{i,s}}{D_{i,s}^{\text{cluster}}}
$$

This is directly analogous to the exon-bin DEXSeq model: the test asks whether the
proportion of site-usage captured by intron $i$ changes between conditions.

---

### 2b. `site_depth` denominator (exon-adjacent window depth)

**Definition**

For each intron $i$ and sample $s$, DSF computes the read depth in a window just
*outside* the intron at each splice-site boundary (the exon-adjacent coverage, not
the intron interior). The denominator is the maximum of the two boundary windows:

$$
D_{i,s}^{\text{site}} = \max\!\left(\text{left\_adjacent\_depth}_{i,s},\; \text{right\_adjacent\_depth}_{i,s}\right)
$$

**The implied $Y_{\text{other}}$**

$$
Y_{\text{other},i,s}^{\text{site}} = D_{i,s}^{\text{site}} - Y_{i,s}
$$

This "other" count represents *reads that were present at the splice-site boundary
but did NOT splice into intron $i$. That pool includes:

- reads that spliced into a different (competing) intron at the same site, and
- reads that read through without splicing (retention / exon-body coverage).

**Important difference from `cluster_max`**: $D_{i,s}^{\text{site}}$ is derived from
continuous coverage depth, not from a sum of discrete junction reads. The "other"
count therefore mixes unspliced and differently-spliced molecules. This makes
$Y_{\text{other}}^{\text{site}}$ *less clean* as an NB-distributed count — it is the
difference of a depth estimate and a junction count, which are measured on different
footprints.

In the DEXSeq paired-count rephrase, the NB assumption on $Y_{\text{other}}^{\text{site}}$
is a weaker approximation than for `cluster_max`. The offset formulation may be
preferred here precisely because it does not require $D^{\text{site}}$ to behave as
a count.

---

### 2c. `splice_plus_retained` denominator

**Definition**

For each splice-site boundary, DSF computes a *splice-plus-retained* depth: the
number of molecules that either spliced at that boundary **or** read continuously
across it into the intron body (i.e., retention-supporting reads). The denominator
is the maximum across the two boundaries:

$$
D_{i,s}^{\text{spr}} = \max\!\left(\text{splice\_depth}_{\text{left}} + \text{retained\_depth}_{\text{left}},\; \text{splice\_depth}_{\text{right}} + \text{retained\_depth}_{\text{right}}\right)
$$

**The implied $Y_{\text{other}}$**

$$
Y_{\text{other},i,s}^{\text{spr}} = D_{i,s}^{\text{spr}} - Y_{i,s}
$$

This "other" count captures two distinct competing signals:

1. reads that spliced at the same boundary into a *different* intron, and
2. reads that were *retained* (not spliced out at this position).

This is the richest denominator: it explicitly acknowledges that a molecule at a
splice site has two choices — splice or retain — and the "other" count represents
all evidence for the alternative outcomes.

In a DEXSeq-style framing this is a natural choice when the biological question is:

> Did the ratio of "spliced via intron $i$" to "not-spliced-at-this-boundary"
> change between conditions?

The paired model then directly tests that competition. The complication is that
the "other" count conflates two biologically distinct events (alternative splicing
and intron retention), so effect-size interpretation requires care.

---

### 2d. `gene_median_splice_plus_retained` denominator

**Definition**

The same splice-plus-retained depth as in 2c, but the per-intron-per-sample value
is replaced by the **median across all introns in the gene** within each sample:

$$
D_{i,s}^{\text{gmed}} = \text{median}_{j \in \text{gene}(i)}\!\left(D_{j,s}^{\text{spr}}\right)
$$

**Motivation**

At very low-count loci, $D_{i,s}^{\text{spr}}$ can be dominated by noise. Sharing
a gene-level median pools information across all the gene's splice sites and
stabilises the denominator. It is similar to how edgeR's TMM normalisation
stabilises library-size estimates by sharing information across genes.

**Implication for the DEXSeq framing**

Because $D_{i,s}^{\text{gmed}}$ is a smoothed estimate (a median of depth values
rather than a directly observed count), $Y_{\text{other},i,s}^{\text{gmed}}$ is not
directly observed — it inherits the smoothing. The paired-count NB model is therefore
an approximation; using it as a fixed offset (the current DSF approach) is actually
more defensible here than constructing a literal "other count" from a smoothed
denominator.

---

### 2e. `splice_plus_retained_betabinom` mode

This mode already moves closest to the DEXSeq paired-count spirit. Rather than
supplying the denominator as an NB offset, it constructs explicit focal and
remainder trials:

$$
\text{successes} = Y_{i,s}, \qquad \text{failures} = D_{i,s}^{\text{spr}} - Y_{i,s}
$$

and fits a **quasi-binomial GLM**. This is equivalent to a DEXSeq-style paired model
under a binomial (rather than NB) distributional assumption. The quasi-binomial
accommodates overdispersion via a variance inflation factor, but it does not model
the mean-variance relationship as flexibly as the NB does.

A full DEXSeq-style rephrase would replace the quasi-binomial with two NB
observations sharing an estimated dispersion — which would be a direct port of the
DEXSeq exon-bin model onto intron features.

---

## 3. Summary comparison

| DSF mode | denominator $D_{i,s}$ | "other" count meaning | DEXSeq rephrase quality |
|---|---|---|---|
| `cluster_max` | max(donor cluster total, acceptor cluster total) | all other junction reads at the same splice site | **Clean**: both focal and other are discrete junction counts; NB assumption is natural |
| `site_depth` | max adjacent exon-window depth | reads at boundary not supporting focal intron (mixed spliced + unspliced) | **Weak**: denominator mixes count types; offset formulation preferred |
| `splice_plus_retained` | max(left spr depth, right spr depth) | reads that spliced elsewhere + reads retained at boundary | **Moderate**: biologically rich; NB on mixture of two event types |
| `gene_median_splice_plus_retained` | per-gene median of spr depth | smoothed "other" — not directly observed | **Approximation**: smoothed denominator breaks the literal paired-count interpretation; keep as offset |
| `betabinom` | max spr depth | explicit focal/rest trials; quasibinomial | **Closest to DEXSeq**: replace quasibinomial with NB for a full DEXSeq-style model |

---

## 4. The singleton-cluster problem and why it matters for the DEXSeq framing

When an intron is the **only** intron using both its donor and its acceptor (a
"both-sided singleton"), the cluster denominator is just the intron itself:

$$
D_{i,s}^{\text{cluster}} = Y_{i,s}
$$

which forces $p_{i,s} = 1$ in every sample, so no differential usage can be
detected. The `cluster_max` code handles this by falling back to `site_depth` for
singletons (`cluster_with_site_singleton_fallback` mode).

In the DEXSeq paired-count framing this manifests as:

$$
Y_{\text{other},i,s}^{\text{cluster}} = D_{i,s}^{\text{cluster}} - Y_{i,s} = 0
$$

A zero "other" count in all samples makes the paired-count model non-identifiable
for that intron. The same fallback strategy (substituting a site-depth or
splice-plus-retained denominator for singletons) is required in the DEXSeq framing
as well.

---

## 5. What a full DEXSeq-style NB model for introns would look like

For each intron $i$ and sample $s$, construct two rows:

$$
Y_{i,s} \sim NB(\mu_{i,s},\; \phi_i)
$$

$$
Y_{\text{other},i,s} \sim NB(\mu_{\text{other},i,s},\; \phi_i)
$$

with a shared dispersion $\phi_i$. The linear predictor for the focal count is:

$$
\log(\mu_{i,s}) = \alpha_i + \beta_i \cdot \text{Group}_s + \delta_s
$$

and symmetrically for the "other" count, where $\delta_s$ is a sample-level
intercept (absorbs sequencing depth). The **condition $\times$ feature interaction**
$\beta_i$ is the test statistic — it measures how the focal intron's share of the
total available evidence shifts between conditions, after accounting for
sample-to-sample variation in that total.

Compared to the offset approach:

| property | offset (current DSF) | DEXSeq paired NB |
|---|---|---|
| $D_{i,s}$ treated as | fixed, known | random; both sides modelled as NB |
| uncertainty in $D_{i,s}$ | ignored | propagated through dispersion |
| information from "other" counts | none | informs dispersion estimation |
| implementation | edgeR with `offset` argument | DEXSeq-style stacked count matrix |
| most suitable $D_{i,s}$ choice | any (depth estimates acceptable) | `cluster_max` (cleanest NB counts); spr modes need care |

The practical difference is largest at **low counts**, where treating $D_{i,s}$ as
fixed underestimates uncertainty and can inflate false positives. At moderate-to-high
counts the two approaches converge.

---

## 6. Relationship to the transcript-level offset model

The offset model described in the companion document
(*GLMs for Differential Transcript Usage with a Gene-Expression Offset*) uses the
parent-gene total $G_i$ as the offset — conceptually the same as `cluster_max` where
the "cluster" is the entire gene. The DEXSeq extension there constructs:

$$
Y_{\text{other},t,i} = G_i - Y_{t,i}
$$

and models both under NB. Every denominator mode in DSF is an intron-level analogue
of that gene-level offset, with varying degrees of locality:

| scope | denominator | analogy |
|---|---|---|
| gene | $G_i = \sum_t Y_{t,i}$ | transcript-level DTU offset |
| splice-site cluster | max(donor total, acceptor total) | intron-level `cluster_max` |
| splice boundary window | adjacent exon depth | `site_depth` |
| splice + retention boundary | spr depth | `splice_plus_retained` |

Moving from gene-level to splice-boundary-level denominators increases specificity
(the "other" becomes a tighter local competitor) at the cost of making the NB
assumption on $Y_{\text{other}}$ less clean.
