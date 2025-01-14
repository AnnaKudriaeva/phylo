# pylint: disable=W0621
import math
import sys

import dash
from dash import Input, Output, html, callback, dcc
import dash_cytoscape as cyto


try:
    from Bio import Phylo
except ModuleNotFoundError as e:
    print(
        e,
        "Please make sure biopython is installed correctly before running this example.",
    )
    sys.exit(1)


def generate_elements(tree, xlen=2, ylen=30, grabbable=False):
    def get_col_positions(tree, column_width=5550):
        """Create a mapping of each clade to its column position."""
        taxa = tree.get_terminals()

        max_label_width = max(len(str(taxon)) for taxon in taxa)
        drawing_width = column_width - max_label_width - 1

        depths = tree.depths()
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)

        fudge_margin = int(math.ceil(math.log(len(taxa), 2)))
        cols_per_branch_unit = (drawing_width - fudge_margin) / float(
            max(depths.values())
        )
        return dict(
            (clade, int(blen * cols_per_branch_unit + 1.0))
            for clade, blen in depths.items()
        )

    def get_row_positions(tree):
        """Calculate row positions for each clade."""
        taxa = tree.get_terminals()
        positions = dict((taxon, 2 * idx) for idx, taxon in enumerate(taxa))

        def calc_row(clade):
            for subclade in clade:
                if subclade not in positions:
                    calc_row(subclade)
            positions[clade] = (
                positions[clade.clades[0]] + positions[clade.clades[-1]]
            ) // 2

        calc_row(tree.root)
        return positions

    def add_to_elements(clade, clade_id):
        """Add nodes and edges to Cytoscape elements."""
        children = clade.clades

        pos_x = col_positions[clade] * xlen
        pos_y = row_positions[clade] * ylen

        cy_source = {
            "data": {"id": clade_id},
            "position": {"x": pos_x, "y": pos_y},
            "classes": "nonterminal",
            "grabbable": grabbable,
        }
        nodes.append(cy_source)

        if clade.is_terminal():
            cy_source["data"]["name"] = clade.name
            cy_source["classes"] = "terminal"

        for n, child in enumerate(children):
            support_id = clade_id + "s" + str(n)
            child_id = clade_id + "c" + str(n)
            pos_y_child = row_positions[child] * ylen

            cy_support_node = {
                "data": {"id": support_id},
                "position": {"x": pos_x, "y": pos_y_child},
                "grabbable": grabbable,
                "classes": "support",
            }

            cy_support_edge = {
                "data": {
                    "source": clade_id,
                    "target": support_id,
                    "sourceCladeId": clade_id,
                },
            }

            cy_edge = {
                "data": {
                    "source": support_id,
                    "target": child_id,
                    "length": clade.branch_length,
                    "sourceCladeId": clade_id,
                },
            }

            if clade.confidence and clade.confidence.value:
                cy_source["data"]["confidence"] = clade.confidence.value

            nodes.append(cy_support_node)
            edges.extend([cy_support_edge, cy_edge])

            add_to_elements(child, child_id)

    col_positions = get_col_positions(tree)
    row_positions = get_row_positions(tree)

    nodes = []
    edges = []

    add_to_elements(tree.clade, "r")

    return nodes, edges


# Define elements, stylesheet, and layout
tree = Phylo.read("my_phylo_2.xml", "phyloxml")
nodes, edges = generate_elements(tree)
elements = nodes + edges

layout = {"name": "preset"}

# List of names to highlight
highlight_names = ["B/PHUKET/3073/2013", "A/flat-faced bat/Peru/033/2010 H18N11", "A/little_yellow-shouldered_bat/Guatemala/164/2009", "A/Jiangxi-Donghu/346/2013", "A/yunnan/0127/2015 H5N6", "A/Australian shelduck/Western Australia/1756/1983 H15N2", "A/chicken/VietNam/NCVD-016/2008 H5N1", "A/Brevig_Mission/1/1918", "A/mallard/Astrakhan/263/1982 H14N15", 
                   "A/Hangzhou/1/2013", "A/duck/Hunan/8-19/2009", "A/pintail_duck/Alberta/114/1979", "A/mallard duck/Alberta/342/1983 H12N1", "A/Missouri/09/2014", "A/chicken/Guangdong/C273/2011 H6N2", "A/black-headed gull/Sweden/5/99 H16N3", "A/black-headed gull/Sweden/1/1999",
                   "A/swine/Hebei/0116/2017","A/Texas/36/1991", "A/New_York/18/2009", "A/thick-billed_murre/Newfoundland/031/2007", "A/Hong_Kong/3239/2008"
]

stylesheet = [
    {
        "selector": ".nonterminal",
        "style": {
            "label": "data(confidence)",
            "background-opacity": 0,
            "text-halign": "left",
            "text-valign": "top",
        },
    },
    {"selector": ".support", "style": {"background-opacity": 0}},
    {
        "selector": "edge",
        "style": {
            "source-endpoint": "inside-to-node",
            "target-endpoint": "inside-to-node",
        },
    },
    {
        "selector": ".terminal",
        "style": {
            "label": "data(name)",
            "width": 10,
            "height": 10,
            "text-valign": "center",
            "text-halign": "right",
            "background-color": "#222222",
        },
    },
]

# Add styles for highlighted names
for name in highlight_names:
    escaped_name = name.replace(" ", r"_")
    stylesheet.append(
        {
            "selector": f'[name *= "{escaped_name}"]',
            "style": {
                "background-color": "red",
                "width": 15,
                "height": 15,
                "border-color": "black",
                "border-width": 2,
            },
        }
    )


# Start the app
app = dash.Dash(__name__, suppress_callback_exceptions=True)
server = app.server

app.layout = html.Div(
    [
        cyto.Cytoscape(
            id="cytoscape",
            elements=elements,
            stylesheet=stylesheet,
            layout=layout,
            style={"height": "95vh", "width": "100%"},
        ),
        html.Button("Save to Folder", id="save-btn"),
        html.Div(id="status-message", style={"margin-top": "10px"}),
    ]
)


@app.callback(
    Output("status-message", "children"),
    Input("save-btn", "n_clicks"),
    prevent_initial_call=True,
)
def save_to_folder(n_clicks):
    # Генерация HTML вручную
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Phylogenetic Tree</title>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.1/cytoscape.min.js"></script>
    </head>
    <body>
        <div id="cy" style="width: 100%; height: 95vh;"></div>
        <script>
            var cy = cytoscape({{
                container: document.getElementById('cy'),
                elements: {elements},
                style: {stylesheet},
                layout: {{ name: 'preset' }}
            }});
        </script>
    </body>
    </html>
    """

    # Сохраняем файл в папку
    filename = "phylogenetic_tree.html"
    with open(filename, "w", encoding="utf-8") as file:
        file.write(html_content)

    # Возвращаем сообщение о статусе
    return f"HTML файл сохранен как: {filename}"


if __name__ == "__main__":
    app.run_server(host="0.0.0.0", port=8080, debug=False)
