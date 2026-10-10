import pytest

PAGE = """# A page

```{{code-cell}} ipython3
{options}x = 1
```

```python
y = x + 1
```
"""


@pytest.mark.parametrize(
    "options",
    ["", ":tags: [hide-input]\n:mystnb-key: value\n", "---\ntags: [hide-input]\nraises-exception: false\n---\n"],
)
def test_run_page_skips_myst_cell_options(run_page, tmp_path, options):
    page = tmp_path / "page.md"
    page.write_text(PAGE.format(options=options), encoding="utf-8")
    assert {name: run_page(page)[name] for name in ("x", "y")} == {"x": 1, "y": 2}


@pytest.mark.parametrize("magic", ["%matplotlib inline", "!pip list"])
def test_run_page_fails_clearly_on_ipython_magics(run_page, tmp_path, magic):
    page = tmp_path / "page.md"
    page.write_text(PAGE.format(options=f"{magic}\n"), encoding="utf-8")
    with pytest.raises(ValueError, match="IPython magic"):
        run_page(page)
