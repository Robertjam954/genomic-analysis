# Re-run after state reset
import pandas as pd
from caas_jupyter_tools import display_dataframe_to_user

path = "/mnt/data/table_1_BBB_regulation.xlsx"
df0 = pd.read_excel(path, header=None)
py
df1 = pd.read_excel(path)

print("Raw shape (no header):", df0.shape)
print("Top 6 rows (no header):")
print(df0.head(6))

print("\nWith default header (row 0 as header):")
print("Columns:", list(df1.columns))
print(df1.head(3))

# Find row that contains 'Layer' as a header cell
row_with_layer = None
for i, row in df0.iterrows():
    if row.astype(str).str.contains(r'(?i)^layer').any():
        row_with_layer = i
        break
print("\nRow index that seems to contain 'Layer' header candidate:", row_with_layer)

if row_with_layer is not None:
    print("\nCandidate header row values:")
    print(df0.iloc[row_with_layer].tolist())
    print("\nRows after header:")
    print(df0.iloc[row_with_layer+1:row_with_layer+6])

# Display preview
display_dataframe_to_user("table_1_BBB_regulation (raw, first 15 rows no-header)", df0.head(15))
