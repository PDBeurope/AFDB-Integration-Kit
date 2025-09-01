"""
PDBFileEditor: A utility for reading, writing, and modifying PDB files.

This module provides a class to handle the fixed-width and multi-line
formatting of PDB file headers and metadata, making it easy to add
and update records like TITLE, REMARK, COMPND, and SOURCE.
"""

import textwrap

class PDBFileEditor:
    """
    A class to read, edit, and write PDB files, with a focus on
    header and metadata records.
    """
    def __init__(self, pdb_content=None):
        """
        Initializes the PDB editor with existing content or an empty list.

        Args:
            pdb_content (list, optional): A list of strings, where each string
                                         is a line of the PDB file. Defaults to None.
        """
        if pdb_content:
            self.lines = pdb_content
        else:
            self.lines = []
        
        # A list to store header/metadata lines to be inserted in order.
        self._header_lines_to_insert = []
        
        # Define the canonical order of PDB records for insertion
        self._record_order = [
            'HEADER', 'TITLE', 'COMPND', 'SOURCE', 'REMARK', 'DBREF',
            'SEQRES', 'CRYST1', 'ORIGX', 'SCALE', 'MODEL', 'ATOM'
        ]

    def load_pdb(self, filename):
        """
        Loads PDB content from a file.

        Args:
            filename (str): The path to the PDB file.
        """
        with open(filename, 'r') as f:
            self.lines = f.readlines()

    def write_pdb(self, filename):
        """
        Writes the current PDB content to a file.
        The method now assembles the final file content by combining the
        header/metadata and the core ATOM lines in the correct order.

        Args:
            filename (str): The path where the PDB file will be saved.
        """
        # Sort the header lines based on their record type to enforce order
        sorted_header_lines = sorted(self._header_lines_to_insert, key=self._get_record_order_key)
        
        # Find the index of the first ATOM record
        atom_start_index = next((i for i, line in enumerate(self.lines) if line.startswith('ATOM')), len(self.lines))
        
        final_lines = []
        # Insert records before ATOM
        for line in sorted_header_lines:
            if self._get_record_name(line) not in ['MODEL']:
                 final_lines.append(line)

        # Insert MODEL record right before ATOM records
        for line in sorted_header_lines:
            if self._get_record_name(line) == 'MODEL':
                final_lines.append(line)
        
        # Append the original ATOM lines
        final_lines.extend(self.lines[atom_start_index:])

        # # Add END record at the end
        # final_lines.append("ENDMDL\n")
        # final_lines.append("END")

        with open(filename, 'w') as f:
            f.writelines(final_lines)
    
    def _get_record_name(self, line):
        """Helper to get the record name from a line."""
        return line.strip().split()[0]
        
    def _get_record_order_key(self, line):
        """Helper to get the sort key for a PDB record line."""
        record_name = self._get_record_name(line)
        try:
            return self._record_order.index(record_name)
        except ValueError:
            return len(self._record_order) # Put unknown records at the end

    def add_title(self, title):
        """
        Adds a TITLE record to the PDB file. Handles multi-line formatting
        and continuation numbers.

        Args:
            title (str): The title of the structure.
        """
        max_width = 70
        wrapped_lines = textwrap.wrap(title, width=max_width)

        for i, line_content in enumerate(wrapped_lines):
            continuation = f"{i + 1: >3}" if i > 0 else "  "
            formatted_line = f"TITLE  {continuation} {line_content:<{max_width}}"
            self._header_lines_to_insert.append(formatted_line + '\n')

    def add_compnd(self, molecule_id, molecule, chain):
        """
        Adds COMPND records for the molecule and chain, supporting multi-line
        molecule descriptions.

        Args:
            molecule_id (int): A unique identifier for the molecule.
            molecule (str): The name of the molecule.
            chain (str): The chain identifier.
        """
        # Line 1: MOL_ID
        self._header_lines_to_insert.append(f"COMPND    MOL_ID: {molecule_id}\n")
        
        # Lines 2 to N: Molecule description (potentially multi-line)
        # Content for MOLECULE starts at column 20 in the AlphaFold example.
        # Line content width is 80 - 19 (COMPND + space + 2-digit number + space + "MOLECULE: ") = 61. 
        molecule_content_prefix = "MOLECULE: "
        wrapped_lines = textwrap.wrap(molecule, width=60)
        
        # First molecule line
        line1 = f"COMPND   2 {molecule_content_prefix}{wrapped_lines[0]};"
        self._header_lines_to_insert.append(f"{line1:<80}\n")
        
        # Subsequent molecule lines, if any
        for i, line_content in enumerate(wrapped_lines[1:]):
            # The continuation number for subsequent lines starts at 3.
            continuation_number = i + 3
            line = f"COMPND  {continuation_number: >2} {line_content};"
            self._header_lines_to_insert.append(f"{line:<80}\n")
        
        # Line N+1: Chain
        # Its continuation number depends on the number of lines the MOLECULE field took.
        chain_continuation_number = 2 + len(wrapped_lines)
        chain_line = f"COMPND  {chain_continuation_number: >2} CHAIN: {chain}"
        self._header_lines_to_insert.append(f"{chain_line:<80}\n")

    def add_source(self, molecule_id, organism_scientific, organism_taxid):
        """
        Adds SOURCE records for the molecule.

        Args:
            molecule_id (int): A unique identifier for the molecule.
            organism_scientific (str): The scientific name of the source organism.
            organism_taxid (int): The NCBI taxonomy ID.
        """
        mol_id_line = f"SOURCE    MOL_ID: {molecule_id};"
        self._header_lines_to_insert.append(f"{mol_id_line:<80}\n")

        organism_line = f"SOURCE   2 ORGANISM_SCIENTIFIC: {organism_scientific};"
        self._header_lines_to_insert.append(f"{organism_line:<80}\n")

        taxid_line = f"SOURCE   3 ORGANISM_TAXID: {organism_taxid}"
        self._header_lines_to_insert.append(f"{taxid_line:<80}\n")

    def add_header(self, pdb_id, date):
        """
        Adds a HEADER record.

        Args:
            pdb_id (str): A PDB identifier.
            date (datetime.date): The date of file creation.
        """
        date_str = date.strftime("%d-%b-%y").upper()
        formatted_line = f"HEADER    {pdb_id:<9}                                  {date_str:<9}                      "
        self._header_lines_to_insert.append(formatted_line + '\n')
    
    def add_model(self, model_number):
        """
        Adds a MODEL record.

        Args:
            model_number (int): The model number.
        """
        formatted_line = f"MODEL     {model_number: >4}"
        self._header_lines_to_insert.append(formatted_line + '\n')
        
    def add_remark_blank_line(self, remark_number=1):
        """
        Adds a blank REMARK line for spacing.
        """
        self._header_lines_to_insert.append(f"REMARK {remark_number: >3} \n")


    def add_remark(self, remark_number, text, numbered=True):
        """
        Adds a REMARK record. Handles general remarks, fixed-width fields,
        and multi-line text wrapping.

        Args:
            remark_number (int): The remark number (e.g., 1 for general remarks).
            text (list or str): A string or list of strings to be added.
            numbered (bool): If True, adds line continuation numbers.
        """
        remark_prefix = f"REMARK {remark_number: >3}"
        
        # If the input is a single string, format it for wrapping.
        if isinstance(text, str):
            text = text.splitlines()

        for line_content in text:
            wrapped_lines = textwrap.wrap(line_content, width=66)
            for i, line_part in enumerate(wrapped_lines):
                continuation_prefix = ""
                if numbered and i > 0:
                    continuation_prefix = f"{i + 1}"
                
                formatted_line = f"{remark_prefix} {continuation_prefix:<2}{line_part:<70}"
                self._header_lines_to_insert.append(formatted_line + '\n')

    def add_remark_reference(self, authors, title, journal, volume, page, year, issn, pmid, doi, reference_number=1):
        """
        Adds a structured REMARK reference section with a specific reference number.

        Args:
            authors (str): String of all authors.
            title (str): Title of the publication.
            journal (str): Journal name.
            volume (str): Journal volume.
            page (str): Starting page number.
            year (str): Publication year.
            issn (str): ISSN.
            pmid (str): PubMed ID.
            doi (str): DOI.
            reference_number (int): The number for this specific reference.
        """
        # Add the 'REFERENCE 1' header line
        self._header_lines_to_insert.append(f"REMARK   1 REFERENCE {reference_number: <2}\n")

        # Format and wrap authors
        author_wrapped = textwrap.wrap(authors, width=61)
        for i, auth_line in enumerate(author_wrapped):
            prefix = f"REMARK   1  AUTH"
            if i > 0:
                prefix += f" {i + 1}"
            formatted_line = f"{prefix:<18} {auth_line:<61}"
            self._header_lines_to_insert.append(formatted_line + '\n')
        
        # Format the title
        title_wrapped = textwrap.wrap(title, width=65)
        for i, titl_line in enumerate(title_wrapped):
            prefix = f"REMARK   1  TITL"
            if i > 0:
                prefix += f" {i + 1}"
            formatted_line = f"{prefix:<14} {titl_line:<65}"
            self._header_lines_to_insert.append(formatted_line + '\n')

        # Format the reference line
        ref_line = f"REMARK   1  REF   {journal:<20} V. {volume: <4} {page: >5} {year: <4}"
        self._header_lines_to_insert.append(f"{ref_line:<80}\n")

        # Format the ISSN line
        issn_line = f"REMARK   1  REFN              ISSN {issn:<9}"
        self._header_lines_to_insert.append(f"{issn_line:<80}\n")
        
        # Format PMID and DOI lines
        pmid_line = f"REMARK   1  PMID  {pmid:<12}"
        self._header_lines_to_insert.append(f"{pmid_line:<80}\n")
        doi_line = f"REMARK   1  DOI   {doi:<66}"
        self._header_lines_to_insert.append(f"{doi_line:<80}\n")


    def add_dbref(self, pdb_id, chain_id, seq_begin, seq_end, db, db_accession, db_id, db_seq_begin, db_seq_end):
        """
        Adds a DBREF record.
        """
        formatted_line = (
            f"DBREF  {pdb_id:<4} {chain_id:<1} {seq_begin: >4} {seq_end: >5} "
            f" {db:<6} {db_accession:<8} {db_id:<12} {db_seq_begin: >5}  {db_seq_end: >5} "
        )
        self._header_lines_to_insert.append(f"{formatted_line:<80}\n")

    def add_crystal_info(self, a, b, c, alpha, beta, gamma):
        """
        Adds CRYST1, ORIGX1, ORIGX2, ORIGX3, SCALE1, SCALE2, SCALE3 records.
        """
        cryst1 = f"CRYST1{a: >9.3f}{b: >9.3f}{c: >9.3f}{alpha: >7.2f}{beta: >7.2f}{gamma: >7.2f} P 1           1"
        self._header_lines_to_insert.append(f"{cryst1:<80}\n")
        
        origx1 = f"ORIGX1{1.000000: >10.6f}{0.000000: >10.6f}{0.000000: >10.6f}{0.00000: >15.5f}"
        self._header_lines_to_insert.append(f"{origx1:<80}\n")
        origx2 = f"ORIGX2{0.000000: >10.6f}{1.000000: >10.6f}{0.000000: >10.6f}{0.00000: >15.5f}"
        self._header_lines_to_insert.append(f"{origx2:<80}\n")
        origx3 = f"ORIGX3{0.000000: >10.6f}{0.000000: >10.6f}{1.000000: >10.6f}{0.00000: >15.5f}"
        self._header_lines_to_insert.append(f"{origx3:<80}\n")
        
        scale1 = f"SCALE1{1.000000: >10.6f}{0.000000: >10.6f}{0.000000: >10.6f}{0.00000: >15.5f}"
        self._header_lines_to_insert.append(f"{scale1:<80}\n")
        scale2 = f"SCALE2{0.000000: >10.6f}{1.000000: >10.6f}{0.000000: >10.6f}{0.00000: >15.5f}"
        self._header_lines_to_insert.append(f"{scale2:<80}\n")
        scale3 = f"SCALE3{0.000000: >10.6f}{0.000000: >10.6f}{1.000000: >10.6f}{0.00000: >15.5f}"
        self._header_lines_to_insert.append(f"{scale3:<80}\n")

    def add_seqres(self, chain_id, sequence_list):
        """
        Adds SEQRES records for a given protein sequence.

        Args:
            chain_id (str): The chain identifier.
            sequence_list (list): A list of three-letter amino acid codes.
        """
        total_residues = len(sequence_list)
        
        # PDB format allows for 13 residues per line.
        num_residues_per_line = 13
        
        for i, start_index in enumerate(range(0, total_residues, num_residues_per_line)):
            serial_num = i + 1
            current_residues = sequence_list[start_index:start_index + num_residues_per_line]
            residue_line_str = " ".join(current_residues)

            # Format the line according to PDB specification
            # Record name (6) + Serial number (3) + Chain ID (1) + Num residues (4) + Residue list (up to 55)
            formatted_line = (
                f"SEQRES {serial_num: >3} {chain_id: >1} {total_residues: >4}  {residue_line_str}"
            )
            self._header_lines_to_insert.append(f"{formatted_line:<80}\n")