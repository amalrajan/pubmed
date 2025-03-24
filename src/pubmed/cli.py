import argparse
from typing import Any, Dict, List, Optional

import pandas as pd
from Bio import Entrez

from pubmed.api import fetch_pubmed_data
from pubmed.parser import parse_article
from pubmed.utils.logger import get_logger


class Args(argparse.Namespace):
    query: str
    file: Optional[str]
    debug: bool
    max: int


def main() -> None:
    """Command-line interface for fetching and analyzing PubMed data.

    Parses arguments, executes search via :func:`fetch_pubmed_data`,
    processes records with :func:`parse_article`, and outputs results
    to CSV file or stdout.

    Usage examples::

        python pubmed_tool.py "cancer[title]" --max 500
        python pubmed_tool.py "machine learning[abstract]" -f output.csv
    """
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description="Fetch PubMed papers with advanced query support",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("query", help="PubMed search query")
    parser.add_argument("-f", "--file", help="Output CSV filename")
    parser.add_argument("-d", "--debug", action="store_true")
    parser.add_argument("--max", type=int, default=100, help="Maximum results to fetch")

    args: Args = parser.parse_args(namespace=Args())
    logger = get_logger(debug=args.debug)

    try:
        logger.debug("Initializing PubMed search with query: %s", args.query)
        records: List[Dict[str, Any]] = fetch_pubmed_data(args.query, args.max)
        logger.info("Successfully retrieved %d records", len(records))

        logger.debug("Parsing article metadata")
        df: pd.DataFrame = pd.DataFrame([parse_article(r) for r in records])

        if args.file:
            logger.info("Writing output to %s", args.file)
            df.to_csv(args.file, index=False)
        else:
            logger.debug("No output file provided. Printing full data to console.")
            print(df.to_csv(index=False))  # Print full CSV to console

    except Entrez.URLError as e:
        logger.error("API Error: %s", str(e))
    except Exception as e:
        logger.critical("Unexpected error: %s", str(e), exc_info=args.debug)
