# import all libraries used in this notebook
import numpy as np
import pandas as pd
import plotnine as p9

def upper_corr_matrix_to_df(draws: np.ndarray) -> pd.DataFrame:
    """
    Compute correlation for a parameter vector.
    Given a 2-d array of draws (rows = draws, columns = variables),
    compute the correlation matrix for the upper triangle.
    - param draws: np.ndarray, N x 1 draws for a single parameter
    - returns a DataFrame ['Var1', 'Var2', 'Correlation']
    """
    cor_mat = np.corrcoef(draws, rowvar=False)

    # Mask the lower triangle (including diagonal)
    upper_cor = cor_mat.copy()
    lower_triangle_indices = np.tril_indices_from(upper_cor)
    upper_cor[lower_triangle_indices] = np.nan

    # Convert to DataFrame and melt
    df_cor = pd.DataFrame(upper_cor)
    melted = df_cor.reset_index().melt(
        id_vars='index',
        var_name='Var2',
        value_name='Correlation'
    )
    melted.rename(columns={'index': 'Var1'}, inplace=True)
    melted = melted.dropna()
    return melted

def plot_icar_corr_matrix(draws: np.ndarray, title: str, size: tuple[int, int]) -> p9.ggplot:
    corr_df = upper_corr_matrix_to_df(draws)
    p = (
        p9.ggplot(corr_df, p9.aes(x='Var1', y='Var2', fill='Correlation'))
        + p9.geom_tile()
        + p9.scale_fill_gradient2(low='darkblue', mid='white', high='darkorange', midpoint=0)
        + p9.theme_minimal()
        + p9.theme(
            figure_size=size,
            axis_text_x=p9.element_blank(),
            axis_text_y=p9.element_blank())
        + p9.ylab('') + p9.xlab('')
        + p9.ggtitle(title)
    )
    return p


def ppc_y_yrep_overlay(
        y_rep: np.ndarray,
        y: pd.Series,
        title: str) -> p9.ggplot:
    y_rep_median = np.median(y_rep, axis=0)
    y_rep_lower = np.percentile(y_rep, 2.5, axis=0)
    y_rep_upper = np.percentile(y_rep, 97.5, axis=0)

    df_plot = pd.DataFrame({
        'obs_id': np.arange(y_rep.shape[1]),
        'y': y,
        'y_rep_median': y_rep_median,
        'y_rep_lower': y_rep_lower,
        'y_rep_upper': y_rep_upper
        })

    # Sort by y, reindex
    df_plot_sorted = df_plot.sort_values(by='y').reset_index(drop=True)
    df_plot_sorted['sorted_index'] = np.arange(len(df_plot_sorted))

    p = (
        p9.ggplot(df_plot_sorted, p9.aes(x='sorted_index'))
        + p9.geom_point(p9.aes(y='y'), color='darkblue', size=0.25)
        + p9.geom_line(p9.aes(y='y_rep_median'), color='orange', alpha=0.8)
        + p9.geom_ribbon(p9.aes(ymin='y_rep_lower', ymax='y_rep_upper'),
                                fill='grey', alpha=0.2)
        + p9.theme(figure_size=(8, 4), axis_text_x=p9.element_blank())
        + p9.ylab('y_rep') + p9.xlab('y')
        + p9.ggtitle(title)
        )
    return(p)


def ppc_central_interval(y_rep: np.ndarray, y: pd.Series) -> str:
    # Compute the 25th and 75th percentiles per observation
    q25 = np.percentile(y_rep, 25, axis=0)
    q75 = np.percentile(y_rep, 75, axis=0)

    # Count how many observed values fall within the central 50% interval
    within_50 = np.sum((y >= q25) & (y <= q75))

    return((
        f"y total: {y_rep.shape[1]}, "
        f"ct y is within y_rep central 50% interval: {within_50}, "
        f"pct: {np.round(100 * within_50 / y_rep.shape[1], decimals=2)}"))


def plot_heatmap(nyc_gdf, data, title, subtitle, scale_name):
    """
    Creates a spatial heatmap
    :param nyc_gdf: GeoDataFrame containing spatial regions
    :param data: Array of values to plot (must match nyc_gdf row order)
    :param title: Plot title
    :param subtitle: Plot subtitle
    :param scale_name: Label for the color scale
    :return: plotnine object
    """
    nyc_gdf = nyc_gdf.copy()
    nyc_gdf["plot_values"] = data
    p = (p9.ggplot(nyc_gdf) +
         p9.geom_map(p9.aes(fill='plot_values')) +
         p9.scale_fill_gradient2(low="blue", mid="white", high="orange", midpoint=0, name=scale_name) +
         p9.labs(title=title, subtitle=subtitle) +
         p9.theme_minimal() +
         p9.theme(figure_size=(20,20),
                  plot_title=p9.element_text(size=32),
                  plot_subtitle=p9.element_text(size=24),
                  legend_position='left',
                  legend_title=p9.element_text(size=20),
                  legend_text=p9.element_text(size=16),
                  legend_key_size=24)
         )
    return p


def ppc_dens_overlay(
    sim_data: pd.Series,
    y_rep: np.ndarray,
    sample_size: int,
    title: str,
    x_label: str
) -> p9.ggplot:
    y_rep_sample = pd.DataFrame(y_rep).sample(n=sample_size).reset_index(drop=True).T
    ppc_dens_plot = p9.ggplot()
    for i in range(sample_size):
        ppc_dens_plot = (ppc_dens_plot
                             + p9.stat_density(mapping=p9.aes(x=y_rep_sample[i]), geom='line', color='lightblue', alpha=0.2))
    ppc_dens_plot = (ppc_dens_plot 
                         + p9.stat_density(mapping=p9.aes(x=sim_data), geom='line', color='darkblue', size=1.1)
                         + p9.ggtitle(title)
                         + p9.xlab(x_label) + p9.ylab("density")
                         + p9.theme(figure_size=(10,5))
         )
    return ppc_dens_plot


def marginal_variances_boxplot(
        names: list[str],
        df_vars: list[pd.DataFrame]
        ) -> p9.ggplot:
    if len(names) != len(df_vars):
        print(f'Error, mismatch between names list and dataframes list')
        return None

    df_list = []
    for model_name, var_array in zip(names, df_vars):
        df_model = pd.DataFrame({
            "index": np.arange(1, df_vars[0].shape[0] + 1),
            "variance": var_array,
            "model": model_name
            })
        df_list.append(df_model)
        df_variances = pd.concat(df_list, ignore_index=True)

    boxplot = (
        p9.ggplot(df_variances, p9.aes(x='model', y='variance', fill='model')) +
        p9.geom_boxplot() +
        p9.ggtitle("Marginal Variances of beta_age Across Models") +
        p9.theme_minimal() +
        p9.labs(x="Model", y="Marginal Variance")
    )
    return boxplot
