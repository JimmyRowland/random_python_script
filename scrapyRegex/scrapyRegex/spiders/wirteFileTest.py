import nbformat

filename = "test2.ipynb"
nb = nbformat.v4.new_notebook()
nb['cells']=[
    # nbformat.v4.new_markdown_cell(response.body),
    nbformat.v4.new_code_cell('''print(1)\nprint(1)\nprint(1)
    '''),
    nbformat.v4.new_code_cell("print(1)")]
nbformat.write(nb,filename)

