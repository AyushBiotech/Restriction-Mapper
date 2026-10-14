# 🧬 Restriction Mapper

A Python-based bioinformatics command-line tool for analyzing DNA sequences using restriction enzymes.

Built using Biopython and Pandas, Restriction Mapper allows you to scan local files or manually pasted sequences, identify restriction cleavage sites, classify enzyme cut patterns, filter results by frequency, and export structured timestamped reports. It also features built-in cross-platform support. 

## ✨ Features

### 1. Flexible DNA Input Methods

- **Manual Entry:** Type or paste your DNA sequence directly into the terminal.
- **File Loading:** Load sequences from local workspace directories with automatic format recognition:
  - `.txt` — Raw DNA sequences
  - `.fasta`, `.fa` — FASTA format
  - `.gb`, `.gbk`, `.genbank` — GenBank format

### 2. Cross-Platform Directory Support

- Tailored startup options for Android, Windows, Linux, and macOS.
- Interactive directory browsing, path validation, and error handling for missing or restricted folders.

### 3. DNA Sequence Validation

Automatically scans input sequences to ensure they contain only standard DNA bases (A, T, G, C), rejecting invalid characters before analysis.

### 4. Topology Selection

Configure restriction analysis for Linear or Circular DNA molecules to match your experimental design.

### 5. Custom or Comprehensive Enzyme Selection

- Test individual restriction enzymes by name.
- Select `All` to evaluate your sequence against the entire catalog via Biopython's `AllEnzymes` collection.
- Build a custom batch of multiple enzymes in a single run.

### 6. Restriction Site Analysis & Classification

For each evaluated enzyme, the program reports:

- Total number of cleavage sites found.
- Exact cutting positions.
- Cut classification types:
  - Blunt ends
  - 5′ overhangs
  - 3′ overhangs
  - Non-cutter

### 7. Cutting-Frequency Filtering

Optionally filter your results by a desired cutting frequency to instantly isolate enzymes that cut your sequence an exact number of times.

### 8. Timestamped CSV Reports

Export your complete analysis results into a structured CSV file named with a precise timestamp (e.g., `restriction_analysis(14-10-2026_19-01-00).csv`).

| Column | Description |
|---|---|
| Enzyme | Restriction enzyme name |
| Total Sites | Number of cleavage sites detected |
| Cutting Positions | Positions of detected cleavage sites |
| Cut Type | Blunt, 5′ overhang, 3′ overhang, or non-cutter |

## 🛠️ Technologies Used

- **Python** — Core programming language
- **Biopython** — Biological sequence handling (`Seq`, `SeqIO`) and restriction enzyme analysis (`RestrictionBatch`, `AllEnzymes`)
- **Pandas** — Data manipulation and structured CSV report generation
- **OS / Datetime modules** — Cross-platform workspace navigation and timestamping

## 📦 Requirements & Installation

Make sure you have Python 3.x installed, then install the required dependencies:

```bash
pip install biopython pandas
```

## 🚀 Getting Started

- Clone or download this repository.
- Run the script in your terminal:

```bash
python restriction_mapper.py
```

- Follow the interactive terminal prompts to select your input source, configure your environment, choose your DNA topology, and select your enzymes.
