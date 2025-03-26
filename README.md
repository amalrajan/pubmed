# PubMed Paper Fetcher

A command-line tool to **fetch and analyze papers from PubMed** using advanced query syntax, and output the results to a CSV file or print to console. The tool supports identifying non-academic authors, affiliated companies, and corresponding author emails.

---

## Features

- **Fetch Papers from PubMed:**  
   Use PubMed's advanced query syntax to retrieve research papers.
- **Flexible Output Options:**  
  - Save results to a CSV file.
  - Print results to the console if no file is specified.
- **Identify Non-Academic Authors & Companies:**  
   Extract and filter authors affiliated with **pharmaceutical/biotech companies**.
- **Error Handling & Debugging:**  
   Enable debug logs to troubleshoot errors during API calls and parsing.

---

## Installation

### Option 1: Install from TestPyPI

```bash
pip install -i https://test.pypi.org/simple/ pubmed-extension
```

### Option 2: Install from Source

```bash
git clone https://github.com/username/pubmed-paper-fetcher.git
cd pubmed-paper-fetcher
poetry install
```

## Usage

```bash
poetry run get-papers-list "cancer AND (biotech[ad] OR pharma[ad])" --debug --file papers.csv
```
