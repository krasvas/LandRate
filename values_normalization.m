%Values normalization function between 0 and 1
function normalized_values=values_normalization(values)
normalized_values=values./max(values);
end