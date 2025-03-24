from typing import Any, Dict, List, Set

from pubmed.utils.logger import get_logger

logger = get_logger()

COMPANY_KEYWORDS: Set[str] = {"pharma", "biotech", "inc", "ltd", "corporation"}
ACADEMIC_KEYWORDS: Set[str] = {"university", "college", "hospital", "institute"}


def parse_article(article: Dict[str, Any]) -> Dict[str, Any]:
    """Parse PubMed MedlineCitation record with industry affiliation detection.

    :param article: Raw article data from Entrez.read()
    :type article: Dict[str, Any]
    :return: Processed article with affiliation analysis
    :rtype: Dict[str, Any]

    Analyzes author affiliations using company indicators (``COMPANY_KEYWORDS``)
    to identify non-academic contributors. Extracts:

    - PubMed ID (PMID)
    - Article title
    - Publication date
    - Semicolon-separated list of non-academic authors
    - Company affiliations
    - Corresponding email if available

    Structure mirrors Biopython's MedlineCitation XML parsing with additional
    industry relationship metadata.
    """
    try:
        medline: Dict[str, Any] = article["MedlineCitation"]
        article_data: Dict[str, Any] = medline.get("Article", {})
        auth_list: List[Dict[str, Any]] = article_data.get("AuthorList", [])

        non_academic_authors: List[str] = []
        company_affiliations: Set[str] = set()

        logger.debug("Processing article with PMID: %s", medline["PMID"])

        for author in auth_list:
            affiliations: List[Dict[str, Any]] = author.get("AffiliationInfo", [])
            for affil in affiliations:
                affil_text: str = affil.get("Affiliation", "").lower()
                if any(kw in affil_text for kw in COMPANY_KEYWORDS):
                    company_affiliations.add(affil_text)
                    non_academic_authors.append(str(author.get("LastName", "")))

        journal_info: Dict[str, Any] = article_data["Journal"]["JournalIssue"]

        parsed_data = {
            "pmid": str(medline["PMID"]),
            "title": str(article_data["ArticleTitle"]),
            "date": journal_info["PubDate"],
            "non_academic_authors": ";".join(non_academic_authors),
            "company_affiliations": ";".join(company_affiliations),
            "corresp_email": str(article_data.get("ELocationID", "")),
        }

        logger.debug("Successfully parsed article with PMID: %s", parsed_data["pmid"])
        return parsed_data

    except KeyError as e:
        logger.error("Missing key in article data: %s", str(e), exc_info=True)
        return {
            "pmid": "N/A",
            "title": "N/A",
            "date": "N/A",
            "non_academic_authors": "",
            "company_affiliations": "",
            "corresp_email": "",
        }
    except Exception as e:
        logger.critical("Error while parsing article: %s", str(e), exc_info=True)
        raise
