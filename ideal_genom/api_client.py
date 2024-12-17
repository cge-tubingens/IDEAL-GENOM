import requests
import time
import logging
import json

class VEPEnsemblRestClient:

    def __init__(self, server='https://rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def _rate_limit(self):
        """Rate-limiting to respect requests per second"""
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(2 - delta)
            self.last_req = time.time()
            self.req_count = 0

    def perform_rest_action(self, method, endpoint, headers=None, params=None, data=None):
        """General method to perform REST actions with GET or POST"""
        self._rate_limit()

        if headers is None:
            headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
        elif 'Content-Type' not in headers:
            headers['Content-Type'] = 'application/json'

        url = self.server + endpoint
        
        try:
            if method == 'GET':
                response = requests.get(url, headers=headers, params=params, timeout=10)
            elif method == 'POST':
                response = requests.post(url, headers=headers, json=data, timeout=10)
            response.raise_for_status()  # Will raise HTTPError for bad responses (4xx, 5xx)
            self.req_count += 1
            return response.json()
        except requests.exceptions.HTTPError as e:
            if response.status_code == 429:
                retry_after = int(response.headers.get('Retry-After', 1))
                logging.warning(f"Rate-limited. Retrying after {retry_after} seconds...")
                time.sleep(retry_after)
                return self.perform_rest_action(method, endpoint, headers, params, data)
            else:
                logging.error(f"HTTP error {response.status_code}: {e}")
        except requests.exceptions.RequestException as e:
            logging.error(f"Request failed: {e}")
        return None

    def post_vep_request(self, ids):
        """Specific method to perform the VEP POST request"""
        
        endpoint = "/vep/human/id"
        headers = {
            "Content-Type": "application/json", 
            "Accept": "application/json"
        }
        data = {
            "ids": ids
        }
        return self.perform_rest_action('POST', endpoint, headers=headers, data=data)

class GeneEnsemblRestClient:
    
    def __init__(self, server='https://rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def _rate_limit(self):
        """Rate-limiting to respect requests per second"""
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

    def perform_rest_action(self, method, endpoint, headers=None, params=None, data=None):
        """General method to perform REST actions with GET or POST"""
        self._rate_limit()

        if headers is None:
            headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
        elif 'Content-Type' not in headers:
            headers['Content-Type'] = 'application/json'

        url = self.server + endpoint
        
        try:
            if method == 'GET':
                response = requests.get(url, headers=headers, params=params, timeout=10)
            elif method == 'POST':
                response = requests.post(url, headers=headers, json=data, timeout=10)
            response.raise_for_status()  # Will raise HTTPError for bad responses (4xx, 5xx)
            self.req_count += 1
            return response.json()
        except requests.exceptions.HTTPError as e:
            if response.status_code == 429:
                retry_after = int(response.headers.get('Retry-After', 1))
                logging.warning(f"Rate-limited. Retrying after {retry_after} seconds...")
                time.sleep(retry_after)
                return self.perform_rest_action(method, endpoint, headers, params, data)
            else:
                logging.error(f"HTTP error {response.status_code}: {e}")
        except requests.exceptions.RequestException as e:
            logging.error(f"Request failed: {e}")
        return None

    def get_gene_location(self, gene_symbol):
        """
        Fetch the start, end, strand, and chromosome for a given gene symbol using the Ensembl REST API.

        Parameters:
            gene_symbol (str): The symbol of the gene (e.g., "BRCA2").

        Returns:
            dict: A dictionary containing start, end, strand, and chromosome information, or an error message.
        """
        # Step 1: Get Ensembl Gene ID from the gene symbol
        ext_xrefs = f"/xrefs/symbol/homo_sapiens/{gene_symbol}"
        xrefs_data = self.perform_rest_action('GET', ext_xrefs)

        if not xrefs_data:
            return {"error": f"Failed to fetch Ensembl ID for gene symbol {gene_symbol}."}

        # Parse response to find the Gene ID
        gene_id = None
        for entry in xrefs_data:
            if entry.get("type") == "gene":
                gene_id = entry.get("id")
                break

        if not gene_id:
            return {"error": f"Gene symbol {gene_symbol} not found in Ensembl."}

        # Step 2: Use the Gene ID to fetch gene details
        ext_lookup = f"/lookup/id/{gene_id}"
        gene_data = self.perform_rest_action('GET', ext_lookup)

        if not gene_data:
            return {"error": f"Failed to fetch gene details for {gene_id}."}

        # Parse and return gene details
        return {
            "gene_symbol": gene_symbol,
            "ensembl_id": gene_id,
            "start": gene_data.get("start"),
            "end": gene_data.get("end"),
            "strand": gene_data.get("strand"),
            "chromosome": gene_data.get("seq_region_name")
        }
