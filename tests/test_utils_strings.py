from tagbrewer.utils import strings

def test_slice_iterator():
    string = "abcdefgh"
    slice_len = 5
    string_slices = [i for i in strings.slice_iterator(string, slice_len)]
    assert string_slices == ["abcde",
                             "bcdef",
                             "cdefg",
                             "defgh"]
