import re
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

        # Exclude papers that only have academic affiliations
        if not company_affiliations:
            logger.debug(
                "Excluding article with PMID %s due to only academic affiliations.",
                medline["PMID"],
            )
            return None

        journal_info: Dict[str, Any] = article_data["Journal"]["JournalIssue"]
        emails_by_author: Dict[str, str] = {}
        for author in article_data.get("AuthorList", []):
            author_name = str(author.get("LastName", ""))
            for aff_info in author.get("AffiliationInfo", []):
                aff_text = aff_info.get("Affiliation", "")
                if "@" in aff_text and "." in aff_text:
                    found_emails = re.findall(
                        r"\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b", aff_text
                    )
                    if found_emails:
                        # If the same author appears in multiple affiliations, join emails
                        if author_name in emails_by_author:
                            # Append additional emails (ensure uniqueness if needed)
                            existing = set(emails_by_author[author_name].split(";"))
                            updated = existing.union(found_emails)
                            emails_by_author[author_name] = ";".join(updated)
                        else:
                            emails_by_author[author_name] = ";".join(found_emails)

        parsed_data = {
            "PubmedID": str(medline["PMID"]),
            "Title": str(article_data["ArticleTitle"]),
            "Publication Date": journal_info["PubDate"],
            "Non-academic Author(s)": ";".join(non_academic_authors),
            "Company Affiliation(s)": ";".join(company_affiliations),
            "Corresponding Author Email": str(emails_by_author),
        }

        logger.debug(
            "Successfully parsed article with PMID: %s", parsed_data["PubmedID"]
        )
        return parsed_data

    except KeyError as e:
        logger.error("Missing key in article data: %s", str(e), exc_info=True)
        return {
            "PubmedID": "N/A",
            "Title": "N/A",
            "Publication Date": "N/A",
            "Non-academic Author(s)": "",
            "Company Affiliation(s)": "",
            "Corresponding Author Email": "",
        }
    except Exception as e:
        logger.critical("Error while parsing article: %s", str(e), exc_info=True)
        raise
