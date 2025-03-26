from typing import Any, Dict, List

from Bio import Entrez

from pubmed.utils.logger import get_logger

logger = get_logger()


def fetch_pubmed_data(query: str) -> List[Dict[str, Any]]:
    """Fetch PubMed records using NCBI's Entrez Programming Utilities.

    :param query: PubMed search query using Entrez syntax
    :type query: str
    :return: List of parsed PubMed records as nested dictionaries
    :rtype: List[Dict[str, Any]]
    :raises Entrez.URLError: For network or API connection issues

    Uses ESearch to get PMIDs with history support, then EFetch to retrieve
    full XML records. Maintains NCBI's rate limiting guidelines through
    Biopython's Entrez module.

    Example::

        records = fetch_pubmed_data("biopython[title]")
    """
    logger.info("Starting PubMed search with query: '%s'", query)

    Entrez.email = "amalrajan.in@gmail.com"

    # ESearch to get PMIDs
    try:
        handle: Any = Entrez.esearch(db="pubmed", term=query, usehistory="y")
        search_results: Dict[str, Any] = Entrez.read(handle)  # type: ignore[attr-defined]
        handle.close()
        logger.debug(
            "ESearch completed successfully, found %s results.",
            search_results.get("Count", 0),
        )
    except Exception as e:
        logger.error("Error during ESearch: %s", str(e), exc_info=True)
        raise

    # EFetch to get full records
    try:
        handle = Entrez.efetch(
            db="pubmed",
            webenv=search_results["WebEnv"],
            query_key=search_results["QueryKey"],
            rettype="xml",
            retmode="xml",
        )
        records: List[Dict[str, Any]] = Entrez.read(handle)["PubmedArticle"]
        handle.close()
        logger.info("Successfully retrieved %d records from PubMed.", len(records))
    except Exception as e:
        logger.error("Error during EFetch: %s", str(e), exc_info=True)
        raise

    return records
