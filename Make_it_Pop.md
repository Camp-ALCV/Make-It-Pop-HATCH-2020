---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.1'
      jupytext_version: 1.2.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

```python
the_books = ["The Book Title", "The Bible", "The Outsiders" , "The River Runs Through It"]
books = []
for book in the_books:
    book = book.replace("The ", "")
    books.append(book)
books
```

```python

```
