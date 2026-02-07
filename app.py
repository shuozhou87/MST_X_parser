from __future__ import annotations

import io
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

try:
    from scipy import stats as scipy_stats
except Exception:  # noqa: BLE001
    scipy_stats = None

CHANNEL_ORDER = ["650", "670", "ratio", "unknown"]


def _channel_sort_key(channel: str) -> int:
    try:
        return CHANNEL_ORDER.index(channel)
    except ValueError:
        return len(CHANNEL_ORDER)


def parse_experiment_and_channel(file_path: Path) -> Tuple[str, str]:
    match = re.match(r"^(?P<exp>.+)-(?P<suffix>650|670|R)$", file_path.stem, flags=re.IGNORECASE)
    if not match:
        return file_path.stem, "unknown"

    experiment = match.group("exp")
    suffix = match.group("suffix").upper()
    channel = {"650": "650", "670": "670", "R": "ratio"}.get(suffix, "unknown")
    return experiment, channel


def infer_core_columns(columns: List[str]) -> Tuple[str, str, str]:
    lower_map = {col: col.strip().lower() for col in columns}

    time_candidates = [col for col, lower in lower_map.items() if "time" in lower]
    capillary_candidates = [
        col
        for col, lower in lower_map.items()
        if "capillary" in lower and ("position" in lower or "pos" in lower or lower == "capillary")
    ]

    if not time_candidates:
        raise ValueError("Cannot find time column.")
    if not capillary_candidates:
        raise ValueError("Cannot find capillary column.")

    time_col = time_candidates[0]
    capillary_col = capillary_candidates[0]

    signal_candidates = [col for col in columns if col not in {time_col, capillary_col}]
    if not signal_candidates:
        raise ValueError("Cannot find signal column.")

    signal_col = signal_candidates[0]
    return time_col, signal_col, capillary_col


def _read_raw_csv_table(file_bytes: bytes) -> pd.DataFrame:
    data = pd.read_csv(io.BytesIO(file_bytes), sep=";")
    if data.shape[1] < 3:
        data = pd.read_csv(io.BytesIO(file_bytes))
    if data.shape[1] < 3:
        raise ValueError(f"Unexpected column count ({data.shape[1]}).")
    return data


def read_raw_csv_from_bytes(file_name: str, file_bytes: bytes) -> pd.DataFrame:
    data = _read_raw_csv_table(file_bytes)
    time_col, signal_col, capillary_col = infer_core_columns(list(data.columns))

    clean = data[[time_col, signal_col, capillary_col]].copy()
    clean.columns = ["time_s", "signal", "capillary"]

    clean["time_s"] = pd.to_numeric(clean["time_s"], errors="coerce")
    clean["signal"] = pd.to_numeric(clean["signal"], errors="coerce")
    clean["capillary"] = pd.to_numeric(clean["capillary"], errors="coerce").round().astype("Int64")
    clean = clean.dropna(subset=["time_s", "signal", "capillary"]).copy()
    clean["capillary"] = clean["capillary"].astype(int)

    experiment, channel = parse_experiment_and_channel(Path(file_name))
    clean["experiment"] = experiment
    clean["channel"] = channel
    clean["source_file"] = file_name
    return clean


def load_raw_dataset_from_uploads(
    uploaded_files: List[object],
) -> Tuple[pd.DataFrame, pd.DataFrame, List[str], Dict[str, bytes]]:
    frames: List[pd.DataFrame] = []
    file_index_rows: List[Dict[str, str]] = []
    bad_files: List[str] = []
    xlsx_payloads: Dict[str, bytes] = {}

    for uploaded in uploaded_files:
        file_name = str(getattr(uploaded, "name", ""))
        payload = uploaded.getvalue()
        suffix = Path(file_name).suffix.lower()

        if suffix in {".xlsx", ".xls"}:
            xlsx_payloads[file_name] = payload
            continue
        if suffix != ".csv":
            continue

        experiment, channel = parse_experiment_and_channel(Path(file_name))
        file_index_rows.append(
            {
                "file": file_name,
                "experiment": experiment,
                "channel": channel,
            }
        )
        try:
            frames.append(read_raw_csv_from_bytes(file_name, payload))
        except Exception as exc:  # noqa: BLE001
            bad_files.append(f"{file_name}: {exc}")

    raw_data = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    file_index = pd.DataFrame(file_index_rows)
    return raw_data, file_index, bad_files, xlsx_payloads


def load_sample_info_from_bytes(file_bytes: bytes) -> pd.DataFrame:
    data = pd.read_excel(io.BytesIO(file_bytes))
    lower_map = {col: col.strip().lower() for col in data.columns}

    capillary_candidates = [col for col, lower in lower_map.items() if "capillary" in lower]
    ligand_candidates = [col for col, lower in lower_map.items() if "ligand" in lower]
    analyte_candidates = [col for col, lower in lower_map.items() if "analyte" in lower]

    if not capillary_candidates:
        raise ValueError("Sample info file does not contain a capillary column.")

    capillary_col = capillary_candidates[0]
    ligand_col = ligand_candidates[0] if ligand_candidates else None
    analyte_col = analyte_candidates[0] if analyte_candidates else None

    out = pd.DataFrame()
    out["capillary"] = pd.to_numeric(data[capillary_col], errors="coerce").round().astype("Int64")
    out["ligand"] = data[ligand_col].astype(str) if ligand_col else "Unknown"
    out["analyte"] = data[analyte_col].astype(str) if analyte_col else "Unknown"

    out = out.dropna(subset=["capillary"]).copy()
    out["capillary"] = out["capillary"].astype(int)
    out["ligand"] = out["ligand"].replace({"nan": "Unknown"})
    out["analyte"] = out["analyte"].replace({"nan": "Unknown"})

    return out.drop_duplicates(subset=["capillary"], keep="first")


def build_sample_group(data: pd.DataFrame, mode: str) -> pd.Series:
    if mode == "Ligand":
        return data["ligand"].fillna("Unknown")
    if mode == "Analyte":
        return data["analyte"].fillna("Unknown")
    if mode == "Capillary position":
        return "Cap-" + data["capillary"].astype(str)
    return data["ligand"].fillna("Unknown") + " | " + data["analyte"].fillna("Unknown")


def add_sample_info(raw_data: pd.DataFrame, sample_info: Optional[pd.DataFrame]) -> pd.DataFrame:
    if sample_info is None or sample_info.empty:
        out = raw_data.copy()
        out["ligand"] = "Unknown"
        out["analyte"] = "Unknown"
        return out

    out = raw_data.merge(sample_info, on="capillary", how="left")
    out["ligand"] = out["ligand"].fillna("Unknown")
    out["analyte"] = out["analyte"].fillna("Unknown")
    return out


def compute_window_signal(data: pd.DataFrame, start_time: float, end_time: float) -> pd.DataFrame:
    window = data[(data["time_s"] >= start_time) & (data["time_s"] <= end_time)].copy()
    if window.empty:
        return pd.DataFrame()

    keys = ["merge_group", "experiment", "channel", "capillary", "sample_group", "ligand", "analyte"]
    per_trace = (
        window.groupby(keys, as_index=False)["signal"]
        .mean()
        .rename(columns={"signal": "window_signal"})
        .sort_values(["channel", "merge_group", "sample_group", "capillary"])
    )
    return per_trace


def pvalue_to_star(p_value: float) -> str:
    if pd.isna(p_value):
        return "NA"
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


def run_reference_ttests(data: pd.DataFrame, reference_group: str) -> pd.DataFrame:
    if scipy_stats is None:
        return pd.DataFrame()

    ref_values = data.loc[data["sample_group"] == reference_group, "window_signal"].dropna()
    if ref_values.empty:
        return pd.DataFrame()

    rows = []
    for group_name in sorted(data["sample_group"].unique().tolist()):
        if group_name == reference_group:
            continue
        group_values = data.loc[data["sample_group"] == group_name, "window_signal"].dropna()
        if group_values.empty:
            continue
        p_value = scipy_stats.ttest_ind(ref_values, group_values, equal_var=False, nan_policy="omit").pvalue
        rows.append(
            {
                "contrast": f"{reference_group} vs {group_name}",
                "group": group_name,
                "n_reference": int(ref_values.shape[0]),
                "n_group": int(group_values.shape[0]),
                "p_value": float(p_value),
                "significance": pvalue_to_star(float(p_value)),
            }
        )
    return pd.DataFrame(rows)


def run_anova(data: pd.DataFrame) -> Optional[float]:
    if scipy_stats is None:
        return None
    group_values = []
    for _, group_df in data.groupby("sample_group"):
        values = group_df["window_signal"].dropna().values
        if values.size > 0:
            group_values.append(values)
    if len(group_values) < 2:
        return None
    return float(scipy_stats.f_oneway(*group_values).pvalue)


def main() -> None:
    st.set_page_config(page_title="Monolith X Raw Data App", layout="wide")
    st.title("NanoTemper Monolith X Raw Data Explorer")
    st.sidebar.markdown("### Upload Files")
    uploaded_files = st.sidebar.file_uploader(
        "Drag CSV/XLSX files here",
        type=["csv", "xlsx"],
        accept_multiple_files=True,
    )
    if not uploaded_files:
        st.info("Upload raw CSV files (for example E2/E3/E4) and optional Sample_Info.xlsx to start.")
        return

    raw_data, file_index, bad_files, xlsx_payloads = load_raw_dataset_from_uploads(uploaded_files)
    if raw_data.empty:
        st.error("No valid CSV raw files were loaded from uploads.")
        return

    xlsx_names = sorted(xlsx_payloads.keys())
    default_xlsx_index = 0
    if "Sample_Info.xlsx" in xlsx_names:
        default_xlsx_index = xlsx_names.index("Sample_Info.xlsx") + 1
    selected_xlsx = st.sidebar.selectbox(
        "Sample info file",
        options=["None"] + xlsx_names,
        index=default_xlsx_index,
    )

    sample_info = None
    if selected_xlsx != "None":
        try:
            sample_info = load_sample_info_from_bytes(xlsx_payloads[selected_xlsx])
        except Exception as exc:  # noqa: BLE001
            st.sidebar.error(f"Sample info load failed: {exc}")

    data = add_sample_info(raw_data, sample_info)

    st.sidebar.markdown("### Filters")
    channels = sorted(data["channel"].dropna().unique().tolist(), key=_channel_sort_key)
    default_channels = [ch for ch in channels if ch != "unknown"] or channels
    selected_channels = st.sidebar.multiselect("Fluorescence channels", channels, default=default_channels)

    experiments = sorted(data["experiment"].dropna().unique().tolist())
    selected_experiments = st.sidebar.multiselect("Experiments", experiments, default=experiments)

    filtered = data[
        data["channel"].isin(selected_channels) & data["experiment"].isin(selected_experiments)
    ].copy()
    if filtered.empty:
        st.warning("No rows left after channel/experiment filters.")
        return

    ligand_options = sorted(filtered["ligand"].dropna().astype(str).unique().tolist())
    selected_ligands = st.sidebar.multiselect("Ligands", ligand_options, default=ligand_options)
    filtered = filtered[filtered["ligand"].isin(selected_ligands)].copy()
    if filtered.empty:
        st.warning("No rows left after ligand filter.")
        return

    analyte_options = sorted(filtered["analyte"].dropna().astype(str).unique().tolist())
    selected_analytes = st.sidebar.multiselect("Analytes", analyte_options, default=analyte_options)
    filtered = filtered[filtered["analyte"].isin(selected_analytes)].copy()
    if filtered.empty:
        st.warning("No rows left after analyte filter.")
        return

    st.sidebar.markdown("### Experiment Merge")
    st.sidebar.caption(
        "Assign E2/E3/E4 to the same merge set to pool their replicates in downstream statistics."
    )
    enable_manual_merge = st.sidebar.checkbox("Enable manual experiment merge", value=False)
    if enable_manual_merge:
        merge_count = st.sidebar.slider("Number of merged sets", min_value=2, max_value=6, value=2)
        merge_set_names = [f"Merge-{i}" for i in range(1, merge_count + 1)]
        exp_to_merge: Dict[str, str] = {}
        for exp in selected_experiments:
            exp_to_merge[exp] = st.sidebar.selectbox(
                f"{exp} ->",
                options=merge_set_names,
                index=0,
                key=f"merge_assign_{exp}",
            )
        filtered["merge_group"] = filtered["experiment"].map(exp_to_merge).fillna("Unassigned")
    else:
        filtered["merge_group"] = filtered["experiment"]

    group_mode = st.sidebar.selectbox(
        "Sample grouping",
        options=["Ligand + Analyte", "Ligand", "Analyte", "Capillary position"],
    )
    filtered["sample_group"] = build_sample_group(filtered, group_mode)

    available_groups = sorted(filtered["sample_group"].dropna().unique().tolist())
    selected_groups = st.sidebar.multiselect("Sample groups to compare", available_groups, default=available_groups)
    if selected_groups:
        filtered = filtered[filtered["sample_group"].isin(selected_groups)].copy()
    else:
        st.warning("No sample groups selected.")
        return

    capillary_options = sorted(filtered["capillary"].dropna().unique().tolist())
    selected_capillaries = st.sidebar.multiselect(
        "Capillary positions",
        capillary_options,
        default=capillary_options,
    )
    filtered = filtered[filtered["capillary"].isin(selected_capillaries)].copy()
    if filtered.empty:
        st.warning("No rows left after sample-group/capillary filters.")
        return

    st.markdown(
        f"Loaded **{len(file_index)}** CSV files with **{len(filtered):,}** points after filters."
    )

    if bad_files:
        with st.expander("Files that failed to parse"):
            for msg in bad_files:
                st.write(msg)

    with st.expander("Discovered raw files"):
        preview = file_index.sort_values(["experiment", "channel", "file"]).reset_index(drop=True)
        st.dataframe(preview, use_container_width=True)

    tab_raw, tab_window = st.tabs(["Raw Time-Signal Series", "Window Single-Point Comparison"])

    with tab_raw:
        st.subheader("Raw Time-Signal Series")
        st.caption("Each line is one capillary trace from one experiment.")

        plot_step = st.slider("Downsample step for plotting", min_value=1, max_value=20, value=1)
        color_by = st.selectbox(
            "Trace color",
            options=["Sample group", "Capillary", "Experiment", "Merge group"],
        )

        color_col = {
            "Sample group": "sample_group",
            "Capillary": "capillary",
            "Experiment": "experiment",
            "Merge group": "merge_group",
        }[color_by]

        plot_data = filtered.sort_values(["experiment", "channel", "capillary", "time_s"]).copy()
        if plot_step > 1:
            plot_data["_row_idx"] = plot_data.groupby(["experiment", "channel", "capillary"]).cumcount()
            plot_data = plot_data[plot_data["_row_idx"] % plot_step == 0].copy()
            plot_data = plot_data.drop(columns=["_row_idx"])

        plot_data["trace_id"] = plot_data["experiment"] + " | Cap-" + plot_data["capillary"].astype(str)
        if color_col == "capillary":
            plot_data["capillary"] = "Cap-" + plot_data["capillary"].astype(str)

        facet_col = "channel" if len(selected_channels) > 1 else None
        fig_raw = px.line(
            plot_data,
            x="time_s",
            y="signal",
            color=color_col,
            line_group="trace_id",
            facet_col=facet_col,
            hover_data=["experiment", "channel", "capillary", "sample_group", "ligand", "analyte"],
        )
        fig_raw.update_layout(height=620)
        st.plotly_chart(fig_raw, use_container_width=True)

        show_group_mean = st.checkbox("Show group mean traces", value=False)
        if show_group_mean:
            mean_data = (
                filtered.groupby(["channel", "time_s", "sample_group"], as_index=False)["signal"]
                .mean()
                .sort_values(["channel", "sample_group", "time_s"])
            )
            fig_mean = px.line(
                mean_data,
                x="time_s",
                y="signal",
                color="sample_group",
                facet_col=("channel" if len(selected_channels) > 1 else None),
            )
            fig_mean.update_layout(height=520)
            st.plotly_chart(fig_mean, use_container_width=True)

    with tab_window:
        st.subheader("Window Single-Point Comparison")
        st.caption("Laser-on window is user defined. Single-point value = mean signal in that window.")

        min_time = float(filtered["time_s"].min())
        max_time = float(filtered["time_s"].max())

        default_start = -4.0 if min_time <= -4.0 <= max_time else min_time
        default_end = -2.0 if min_time <= -2.0 <= max_time and -2.0 > default_start else min(default_start + 1.5, max_time)
        if default_end <= default_start:
            default_end = min(max_time, default_start + 0.1)

        window_start, window_end = st.slider(
            "Laser on window (seconds)",
            min_value=min_time,
            max_value=max_time,
            value=(float(default_start), float(default_end)),
            step=0.01,
        )

        per_trace = compute_window_signal(filtered, window_start, window_end)
        if per_trace.empty:
            st.warning("No points inside selected window.")
            return

        col1, col2, col3 = st.columns(3)
        col1.metric("Window start (s)", f"{window_start:.2f}")
        col2.metric("Window end (s)", f"{window_end:.2f}")
        col3.metric("Single-point traces", f"{len(per_trace):,}")

        compare_groups = sorted(per_trace["sample_group"].dropna().unique().tolist())
        default_a = compare_groups[:1]
        default_b = compare_groups[1:2] if len(compare_groups) > 1 else []

        col_a, col_b = st.columns(2)
        group_a = col_a.multiselect("Comparison set A", compare_groups, default=default_a)
        group_b = col_b.multiselect("Comparison set B", compare_groups, default=default_b)

        per_trace["compare_set"] = "Other"
        per_trace.loc[per_trace["sample_group"].isin(group_a), "compare_set"] = "A"
        per_trace.loc[per_trace["sample_group"].isin(group_b), "compare_set"] = "B"

        compare_data = per_trace[per_trace["compare_set"].isin(["A", "B"])].copy()
        if compare_data.empty:
            st.info("Select at least one sample group in set A or set B.")
        else:
            fig_cmp = px.box(
                compare_data,
                x="compare_set",
                y="window_signal",
                color="merge_group",
                facet_col=("channel" if len(selected_channels) > 1 else None),
                points="all",
                hover_data=["experiment", "capillary", "sample_group", "ligand", "analyte"],
            )
            fig_cmp.update_layout(height=560)
            st.plotly_chart(fig_cmp, use_container_width=True)

            compare_stats = (
                compare_data.groupby(["channel", "merge_group", "compare_set"])["window_signal"]
                .agg(count="count", mean="mean", std="std", median="median")
                .reset_index()
            )
            st.dataframe(compare_stats, use_container_width=True)

        st.markdown("### Bar Chart + Significance")
        stats_channels = sorted(per_trace["channel"].dropna().unique().tolist(), key=_channel_sort_key)
        default_channel_index = stats_channels.index("ratio") if "ratio" in stats_channels else 0
        stats_channel = st.selectbox(
            "Channel for statistics",
            options=stats_channels,
            index=default_channel_index,
        )
        stats_input = per_trace[per_trace["channel"] == stats_channel].copy()

        merge_options = sorted(stats_input["merge_group"].dropna().unique().tolist())
        selected_merge_sets = st.multiselect(
            "Merged sets included in bar/stats",
            merge_options,
            default=merge_options,
        )
        stats_input = stats_input[stats_input["merge_group"].isin(selected_merge_sets)].copy()

        if stats_input.empty:
            st.warning("No single-point traces left for selected channel/merge sets.")
        else:
            summary = (
                stats_input.groupby("sample_group")["window_signal"]
                .agg(n="count", mean="mean", std="std")
                .reset_index()
                .sort_values("sample_group")
            )
            summary["sem"] = summary["std"] / summary["n"].pow(0.5)
            error_mode = st.selectbox("Error bar", options=["SEM", "STD"], index=0)
            summary["error"] = summary["sem"] if error_mode == "SEM" else summary["std"]

            group_names = summary["sample_group"].tolist()
            default_ref = 0
            for idx, name in enumerate(group_names):
                if "blank" in str(name).lower():
                    default_ref = idx
                    break
            reference_group = st.selectbox(
                "Reference group for pairwise t-test",
                options=group_names,
                index=default_ref,
            )

            anova_p = run_anova(stats_input)
            ttest_table = run_reference_ttests(stats_input, reference_group)
            if not ttest_table.empty:
                ttest_table["p_value"] = ttest_table["p_value"].map(lambda x: f"{x:.3e}")
            star_map = {reference_group: "Ref"}
            if not ttest_table.empty:
                star_map.update(dict(zip(ttest_table["group"], ttest_table["significance"])))
            summary["sig"] = summary["sample_group"].map(star_map).fillna("")

            fig_bar = go.Figure()
            fig_bar.add_trace(
                go.Bar(
                    x=summary["sample_group"],
                    y=summary["mean"],
                    error_y={"type": "data", "array": summary["error"], "visible": True},
                    text=summary["sig"],
                    textposition="outside",
                    name="Window mean",
                    marker_color="#4C78A8",
                    opacity=0.88,
                )
            )
            fig_bar.add_trace(
                go.Scatter(
                    x=stats_input["sample_group"],
                    y=stats_input["window_signal"],
                    mode="markers",
                    name="Replicates",
                    marker={"size": 8, "color": "rgba(30,30,30,0.55)"},
                )
            )
            fig_bar.update_layout(
                xaxis_title="Sample group",
                yaxis_title=f"Window signal ({stats_channel})",
                height=560,
            )
            st.plotly_chart(fig_bar, use_container_width=True)

            if scipy_stats is None:
                st.warning("`scipy` is not installed, so ANOVA/t-test is unavailable.")
            else:
                if anova_p is not None:
                    st.write(f"One-way ANOVA p-value: `{anova_p:.3e}` ({pvalue_to_star(anova_p)})")
                if ttest_table.empty:
                    st.info("Pairwise t-test result is empty (need at least two groups).")
                else:
                    st.dataframe(ttest_table, use_container_width=True)

        st.markdown("**All single-point results (window mean by trace)**")
        st.dataframe(per_trace, use_container_width=True)

        download_bytes = per_trace.to_csv(index=False).encode("utf-8")
        st.download_button(
            label="Download window single-point CSV",
            data=download_bytes,
            file_name="window_single_point_results.csv",
            mime="text/csv",
        )


if __name__ == "__main__":
    main()
