from bs4 import BeautifulSoup as BS

def parse_allele_frequencies(path):
  """Parses an allele frequency html from the allele frequency database."""
  with open(path, "rb") as file:
    soup = BS(file, "html5lib")
    table = soup.find(class_="tblNormal").contents[1]
    # print(table.contents)
    table_rows = []
    for elem in table:
      if not isinstance(elem, str):
        table_rows.append(elem)
    headers = []
    rows = []
    for child in table_rows[0]:
      if not isinstance(child, str):
        headers.append(" ".join(child.text.replace("\n", " ").split()))
    for idx in range(1, len(table_rows)):
      table_row = table_rows[idx]
      new_row = []
      for child in table_row:
        if not isinstance(child, str):
          new_row.append(" ".join(child.text.replace("\n", " ").replace(",", "").split()))
      rows.append(new_row)
    rows.insert(0, headers)
  return rows

def load_and_dump(path, output):
  """Loads an allele frequency html file and dumps its contents to a csv file."""
  rows = parse_allele_frequencies(path)
  with open(output, "w") as csv_file:
    for row in rows:
      csv_file.write(";".join(row) + "\n")

if __name__ == "__main__":
  import sys
  load_and_dump(sys.argv[1], sys.argv[1] + ".csv")
