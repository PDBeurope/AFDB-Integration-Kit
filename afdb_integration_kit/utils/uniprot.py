import requests
from typing import Dict, Any
import logging

logger = logging.getLogger("afdb_integration_kit")

UNIPROT_API_BASE_URL = "https://rest.uniprot.org/uniprotkb"

class UniprotAPIClient:
    """Client for fetching data from the UniProt API."""

    def __init__(self, session: requests.Session, base_url: str = UNIPROT_API_BASE_URL):
        self.session = session
        self.base_url = base_url

    def fetch_metadata(self, uniprot_accession: str) -> Dict[str, Any]:
        """Fetches metadata for a given UniProt ID."""
        url = f"{self.base_url}/{uniprot_accession}.json"
        logger.info(f"Fetching UniProt metadata for ID: {uniprot_accession}")
        try:
            response = self.session.get(url, timeout=15)
            response.raise_for_status()
            return response.json()
        except requests.RequestException as e:
            logger.warning(
                f"Failed to fetch UniProt ID {uniprot_accession}: {e}. "
                "Please check the ID and your network connection."
            )
            return {}