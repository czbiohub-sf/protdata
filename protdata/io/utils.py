def cleanup_obsvar(df):
    """
    Convert object columns to strings to avoid h5ad serialization issues.
    """
    for col in df.columns:
        if df[col].dtype == "object":
            df[col] = df[col].astype(str)
    return df
