#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import pandas as pd
import plotly.express as px
from dash import Dash, Input, Output, State, dash_table, dcc, html


def build_app(df: pd.DataFrame, export_dir: Path) -> Dash:
    app = Dash(__name__)
    numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    default_x = "rfd3_score" if "rfd3_score" in df.columns else (numeric_cols[0] if numeric_cols else "design_id")
    default_y = "fbc_pass1_score" if "fbc_pass1_score" in df.columns else (numeric_cols[1] if len(numeric_cols) > 1 else default_x)
    hover_cols = [c for c in ["design_id", "charge_net", "fbc_revalidation_pass"] if c in df.columns]

    app.layout = html.Div(
        [
            html.H2("Binder Candidate Selection"),
            html.Div(
                [
                    html.Label("X metric"),
                    dcc.Dropdown(options=df.columns.tolist(), value=default_x, id="x-metric"),
                    html.Label("Y metric"),
                    dcc.Dropdown(options=df.columns.tolist(), value=default_y, id="y-metric"),
                    html.Label("Net charge range"),
                    dcc.RangeSlider(
                        id="charge-range",
                        min=float(df["charge_net"].min()) if "charge_net" in df.columns else -50.0,
                        max=float(df["charge_net"].max()) if "charge_net" in df.columns else 50.0,
                        value=[
                            float(df["charge_net"].min()) if "charge_net" in df.columns else -10.0,
                            float(df["charge_net"].max()) if "charge_net" in df.columns else 10.0,
                        ],
                        step=0.5,
                    ),
                ],
                style={"maxWidth": "900px"},
            ),
            dcc.Graph(id="scatter"),
            html.Button("Export selected", id="export-btn", n_clicks=0),
            html.Div(id="export-status"),
            dash_table.DataTable(
                id="table",
                columns=[{"name": c, "id": c} for c in df.columns],
                data=df.to_dict("records"),
                filter_action="native",
                sort_action="native",
                row_selectable="multi",
                selected_rows=[],
                page_size=20,
            ),
        ]
    )

    @app.callback(
        Output("scatter", "figure"),
        Input("x-metric", "value"),
        Input("y-metric", "value"),
        Input("charge-range", "value"),
    )
    def update_plot(x_metric, y_metric, charge_range):
        dff = df
        if "charge_net" in dff.columns:
            lo, hi = charge_range
            dff = dff[(dff["charge_net"] >= lo) & (dff["charge_net"] <= hi)]
        return px.scatter(
            dff,
            x=x_metric,
            y=y_metric,
            hover_data=hover_cols,
            color="fbc_revalidation_pass" if "fbc_revalidation_pass" in dff.columns else None,
        )

    @app.callback(
        Output("export-status", "children"),
        Input("export-btn", "n_clicks"),
        State("table", "derived_virtual_data"),
        State("table", "derived_virtual_selected_rows"),
        prevent_initial_call=True,
    )
    def export_selected(n_clicks, rows, selected_rows):
        if not rows or not selected_rows:
            return "No rows selected."
        selected = [rows[i] for i in selected_rows]
        export_dir.mkdir(parents=True, exist_ok=True)
        csv_path = export_dir / "shortlist.csv"
        json_path = export_dir / "shortlist.json"
        pd.DataFrame(selected).to_csv(csv_path, index=False)
        json_path.write_text(json.dumps(selected, indent=2))
        return f"Exported {len(selected)} candidates to {csv_path} and {json_path}"

    return app


def main() -> None:
    parser = argparse.ArgumentParser(description="Plotly dashboard for candidate selection.")
    parser.add_argument("--aggregate-csv", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8050)
    args = parser.parse_args()

    df = pd.read_csv(args.aggregate_csv)
    app = build_app(df, args.outdir)
    app.run(host=args.host, port=args.port, debug=False)


if __name__ == "__main__":
    main()

