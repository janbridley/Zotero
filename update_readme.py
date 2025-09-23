import sqlite3
import os
import re
import argparse
from typing import List, Optional
from contextlib import closing

HEADER_TEXT = """
A database of my zotero library, which can be updated by the script in
[update_readme.py](update_readme.py). This allows for lightweight git backups of the
content of a zotero library without burning LFS storage or encountering any copyright
issues. 
"""

SQL_ITEM_METADATA = """
SELECT f.fieldName, idv.value
FROM itemData id
JOIN itemDataValues idv ON id.valueID = idv.valueID
JOIN fields f ON id.fieldID = f.fieldID
WHERE id.itemID = ? AND f.fieldName IN ('title', 'date', 'DOI')
"""
SQL_ITEM_AUTHORS = """
SELECT c.lastName
FROM itemCreators ic
JOIN creators c ON ic.creatorID = c.creatorID
JOIN creatorTypes ct ON ic.creatorTypeID = ct.creatorTypeID
WHERE ic.itemID = ? AND ct.creatorType = 'author'
ORDER BY ic.orderIndex
"""
SQL_ITEM_ATTACHMENT = "SELECT path FROM itemAttachments WHERE parentItemID = ? AND linkMode = 2"
SQL_SUBCOLLECTIONS = "SELECT collectionID, collectionName FROM collections WHERE parentCollectionID = ?"
SQL_SUBCOLLECTIONS_ROOT = "SELECT collectionID, collectionName FROM collections WHERE parentCollectionID IS NULL"
SQL_COLLECTION_ITEMS = """
SELECT i.itemID FROM collectionItems ci
JOIN items i ON ci.itemID = i.itemID
WHERE ci.collectionID = ? AND i.itemTypeID NOT IN (1, 14)
"""
SQL_UNFILED_ITEMS = """
SELECT i.itemID FROM items i
LEFT JOIN collectionItems ci ON i.itemID = ci.itemID
WHERE ci.collectionID IS NULL AND i.itemTypeID NOT IN (1, 14)
"""

def get_year_from_date(date_str: Optional[str]) -> str:
    """Extract the four-digit year from a date string."""
    if not date_str:
        return ""
    match = re.search(r'\d{4}', date_str)
    return match.group(0) if match else ""

def get_authors(cursor: sqlite3.Cursor, item_id: int) -> str:
    """Retrieve and format a comma-separated list of author last names."""
    cursor.execute(SQL_ITEM_AUTHORS, (item_id,))
    authors = cursor.fetchall()
    if not authors:
        return ""
    last_names = [row[0] for row in authors if row[0]]
    return ", ".join(last_names)

def get_item_markdown(cursor: sqlite3.Cursor, item_id: int, storage_folder_name: str) -> Optional[str]:
    """Fetch metadata and format it as a Markdown line, returning None for invalid items."""
    cursor.execute(SQL_ITEM_METADATA, (item_id,))
    metadata = dict(cursor.fetchall())
    
    title = metadata.get('title')
    if not title:
        return None

    cleaned_title = title.strip()
    if cleaned_title in ["PDF", "Full Text"]:
        return None

    cleaned_title = re.sub(r'\.pdf$', '', cleaned_title, flags=re.IGNORECASE)

    authors_str = get_authors(cursor, item_id)
    year = get_year_from_date(metadata.get('date', ''))

    title_part = f'"{cleaned_title}"'
    main_parts = [p for p in [title_part, authors_str] if p]
    main_text = ", ".join(main_parts)
    link_text = f"{main_text} ({year})" if year else main_text

    doi = metadata.get('DOI')
    if doi and doi.strip():
        doi_url = f"https://doi.org/{doi.strip()}"
        return f"[{link_text}]({doi_url})"

    cursor.execute(SQL_ITEM_ATTACHMENT, (item_id,))
    attachment = cursor.fetchone()
    if attachment and attachment[0] and attachment[0].startswith('storage:'):
        relative_path = attachment[0].replace('storage:', f'{storage_folder_name}/')
        return f"[{link_text}]({relative_path})"

    return link_text

def process_collections_recursive(
    cursor: sqlite3.Cursor,
    parent_id: Optional[int],
    level: int,
    storage_folder_name: str
) -> List[str]:
    """Recursively process collections and their items to build Markdown headers."""
    markdown_output: List[str] = []
    header_level = level + 2  # Start with H2 for top-level collections

    query, params = (SQL_SUBCOLLECTIONS_ROOT, ()) if parent_id is None else (SQL_SUBCOLLECTIONS, (parent_id,))
    cursor.execute(query, params)
    collections = cursor.fetchall()

    for collection_id, collection_name in sorted(collections, key=lambda x: x[1]):
        if markdown_output or level > 0:
            markdown_output.append("") # Add space between sections

        markdown_output.append(f"{'#' * header_level} {collection_name}")

        cursor.execute(SQL_COLLECTION_ITEMS, (collection_id,))
        item_lines = []
        for (item_id,) in cursor.fetchall():
            markdown_line = get_item_markdown(cursor, item_id, storage_folder_name)
            if markdown_line:
                item_lines.append(f"- {markdown_line}")
        if item_lines:
            markdown_output.extend(item_lines)

        markdown_output.extend(
            process_collections_recursive(cursor, collection_id, level + 1, storage_folder_name)
        )
    return markdown_output

def generate_readme(zotero_dir: str, readme_name: str, dry_run: bool) -> None:
    """Generate a README file for a Zotero library."""
    zotero_dir = os.path.abspath(zotero_dir)
    if not os.path.isdir(zotero_dir):
        raise ValueError(f"Invalid directory: {zotero_dir}")

    db_path = os.path.join(zotero_dir, 'zotero.sqlite')
    if not os.path.isfile(db_path):
        raise FileNotFoundError(f"Database not found at {db_path}")

    storage_folder_name = 'storage'
    
    with closing(sqlite3.connect(f'file:{db_path}?mode=ro', uri=True)) as conn:
        cursor = conn.cursor()
        markdown_lines = process_collections_recursive(cursor, None, 0, storage_folder_name)

        cursor.execute(SQL_UNFILED_ITEMS)
        unfiled_item_ids = [row[0] for row in cursor.fetchall()]
        if unfiled_item_ids:
            markdown_lines.append("")
            markdown_lines.append("## Unfiled Items")
            for item_id in unfiled_item_ids:
                markdown_line = get_item_markdown(cursor, item_id, storage_folder_name)
                if markdown_line:
                    markdown_lines.append(f"- {markdown_line}")

    output_path = os.path.join(zotero_dir, readme_name)
    output_content = "# Zotero Library\n\n" + HEADER_TEXT + '\n'.join(markdown_lines)
    
    if dry_run:
        print(output_content)
        return

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(output_content)
    print(f"Successfully updated {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a README.md for a Zotero library.")
    parser.add_argument("zotero_directory", help="Path to your Zotero data directory.")
    parser.add_argument("--output", default="README.md", help="Name of the output file (default: README.md).")
    parser.add_argument("--dry-run", action="store_true", help="Preview output without writing a file.")
    args = parser.parse_args()

    try:
        generate_readme(args.zotero_directory, args.output, args.dry_run)
    except Exception as e:
        print(f"Error: {e}")
