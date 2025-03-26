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


def main() -> None:
    """Command-line interface for fetching and analyzing PubMed data.

    Parses arguments, executes search via :func:`fetch_pubmed_data`,
    processes records with :func:`parse_article`, and outputs results
    to CSV file or stdout.

    Usage examples::

        python pubmed_tool.py "cancer[title]
        python pubmed_tool.py "machine learning[abstract]" -f output.csv
    """
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description="Fetch PubMed papers with advanced query support",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("query", help="PubMed search query")
    parser.add_argument("-f", "--file", help="Output CSV filename")
    parser.add_argument("-d", "--debug", action="store_true")

    args: Args = parser.parse_args(namespace=Args())
    logger = get_logger(debug=args.debug)

    try:
        logger.debug("Initializing PubMed search with query: %s", args.query)
        records: List[Dict[str, Any]] = fetch_pubmed_data(args.query)
        logger.info("Successfully retrieved %d records", len(records))

        logger.debug("Parsing article metadata")

        # Skip records with only academic affiliations
        df = pd.DataFrame(
            [
                result
                for record in records
                if (result := parse_article(record)) is not None
            ]
        )

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
