# LaTeX / Math Formatting Guide for GitHub-Compatible Markdown

This guide is written for an AI agent authoring or editing Markdown files that contain
mathematical notation intended to render correctly on **GitHub** (and secondarily in
VS Code's built-in Markdown preview).

---

## 1. Delimiters

GitHub uses its own math extension on top of CommonMark/GFM. Only two delimiter styles
are recognized.

| Style | Delimiter | Renders as |
|---|---|---|
| Inline math | `$...$` | Inline equation |
| Display (block) math | `$$` on its own line, content, then `$$` on its own line | Centered block |

### Do NOT use

- `\(...\)` — not recognized by GitHub; renders as literal text.
- `\[...\]` — not recognized by GitHub; renders as literal text.
- `` `$...$` `` — backticks suppress math; renders as code.

### Correct block math format

```
$$
e^{i\pi} + 1 = 0
$$
```

The opening `$$` and closing `$$` must each be **alone on their own line** with no
other characters (not even trailing spaces).

---

## 2. The Setext Heading Trap (most common GitHub math breakage)

GitHub's Markdown parser runs **before** the math extension. CommonMark defines a
**setext heading**: any line of text followed immediately by a line of `=` or `-`
characters becomes an `<h1>` or `<h2>` heading. This destroys any math block that
contains a lone `=` or `-` operator line.

### Broken (renders as headings, not math on GitHub)

```
$$
\text{usage}_{T1}
=
\frac{80}{100}
=
0.8
$$
```

GitHub sees the `\text{usage}_{T1}` line followed by a lone `=` and promotes it to an
`<h1>`, shattering the block before the math extension ever runs.

### Fixed (keep `=` on the same line as an operand)

```
$$
\text{usage}_{T1} = \frac{80}{100} = 0.8
$$
```

### Rule

> Inside a `$$` block, **never place `=` or `-` alone on its own line**.
> Always attach the operator to an expression on the same line.

A lone `-` is equally dangerous because it also triggers a setext heading (`<h2>`).

---

## 3. Blank Lines Inside Math Blocks

A blank line inside a `$$` block terminates the block in GitHub's parser.

### Broken

```
$$
a = b

c = d
$$
```

### Fixed

```
$$
a = b \\
c = d
$$
```

or split into two separate blocks.

---

## 4. Inline Math in Tables and Blockquotes

GitHub recognizes inline `$...$` inside table cells and blockquotes. Use it freely:

```markdown
| term | meaning |
|---|---|
| $\alpha_t$ | baseline log usage |
| $\beta_t$ | log fold change in usage |
```

```markdown
> The model estimates $\beta_t \approx -0.693$.
```

---

## 5. VS Code vs GitHub Differences

| Feature | VS Code preview | GitHub |
|---|---|---|
| `$...$` inline | ✅ | ✅ |
| `$$...$$` block | ✅ | ✅ |
| `\(...\)` inline | ✅ (KaTeX) | ❌ renders as text |
| `\[...\]` block | ✅ (KaTeX) | ❌ renders as text |
| Lone `=` inside `$$` | ✅ (KaTeX ignores setext) | ❌ breaks block |
| Blank line inside `$$` | ✅ (KaTeX ignores) | ❌ terminates block |

This means a file can look perfectly formatted in VS Code while being broken on GitHub.
**Always verify against GitHub's parser when publishing.**

---

## 6. How to Verify Programmatically

Use GitHub's public Markdown API to render the file and check whether all math spans
are recognized. GitHub wraps recognized math in `<math-renderer>` tags. Any remaining
`$` signs in the output indicate unrecognized (broken) math.

```python
import json, re, subprocess

src = open('your_file.md').read()
payload = json.dumps({'text': src, 'mode': 'gfm'})
open('/tmp/payload.json', 'w').write(payload)

html = subprocess.run(
    ['curl', '-s', '-X', 'POST',
     '-H', 'Content-Type: application/json',
     '--data', '@/tmp/payload.json',
     'https://api.github.com/markdown'],
    capture_output=True, text=True
).stdout

# Strip recognized math
stripped = re.sub(r'<math-renderer.*?</math-renderer>', '', html, flags=re.S)

# Flag any remaining $
for line in stripped.split('\n'):
    if '$' in line:
        print('UNRECOGNIZED:', line.strip()[:200])
```

Run this before every commit that changes math-heavy content.

---

## 7. Quick Reference Checklist

Before committing a math-heavy Markdown file to GitHub:

- [ ] All inline math uses `$...$` (not `\(...\)`)
- [ ] All block math uses `$$` alone on its own line (open and close)
- [ ] No lone `=` or `-` operator lines inside `$$` blocks
- [ ] No blank lines inside `$$` blocks
- [ ] Verified with the GitHub Markdown API script above

---

## 8. Common LaTeX Commands That Work on Both Platforms

These are all safe to use:

```
\frac{a}{b}          fractions
\left( ... \right)   auto-sized parentheses
\begin{cases}        piecewise expressions
\text{...}           upright text inside math
\log \exp \sum \prod standard operators
\alpha \beta \mu     Greek letters
\rightarrow \approx  arrows and relations
\underbrace{...}_{} labelled braces
\times \cdot         multiplication
e^{...}              exponentiation
\qquad \quad         spacing (safe when not alone on a line)
```

The only spacing commands that become dangerous are those that, when combined with
surrounding blank lines, leave only `=` or `-` exposed on a line.
