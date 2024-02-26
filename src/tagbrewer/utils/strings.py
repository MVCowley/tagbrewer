def slice_iterator(lst, slice_len):
    """
    Source: https://stackoverflow.com/questions/1335392/iteration-over-list-slices
    """
    for i in range(len(lst) - slice_len + 1):
        yield lst[i:i + slice_len]