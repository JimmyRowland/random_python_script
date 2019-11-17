import scrapy
import nbformat
from selenium import webdriver


def get_forms(response):
    # return response.css("form")
    return response.xpath("//table")


def get_excise_name(form):
    return form.css("td::text").get()


def get_regex(form):
    return form.xpath("//input/@value").get()


def get_test_string_names(form):
    return form.xpath("//textarea/text()").get().split("\n")


class Lab7Spider(scrapy.Spider):
    name = "regexPage"

    def __init__(self):
        self.driver = webdriver.Firefox()

    def start_requests(self):
        urls = ["file://www.students.cs.ubc.ca/~cs-121/2019W1/Labs/Lab7/lab7.html"]
        for url in urls:
            yield scrapy.Request(url=url, callback=self.parse)

    def parse(self, response):
        page = response.url.split("/")[-4]
        filename = "lab7_regex.ipynb"
        nb = nbformat.v4.new_notebook()
        nb['cells'] = []
        print("!!!!!!!!!!")
        print(get_forms(response))
        print(response.css("table"))

        for form in get_forms(response):
            # nb['cells'].append(nbformat.v4.new_markdown_cell(form.getall()))

            # print(get_excise_name(form))
            nb['cells'].append(nbformat.v4.new_markdown_cell("#" + get_excise_name(form)))
            # print(get_regex(form))
            # nb['cells'].append(nbformat.v4.new_markdown_cell(get_regex(form)))
            # print(get_test_string_names(form))
            # nb['cells'].append(nbformat.v4.new_markdown_cell("\n".join(get_test_string_names(form))))

            nb['cells'].append(nbformat.v4.new_code_cell('''print(1)\nprint(1)'''))
        nb['cells'].append(nbformat.v4.new_code_cell('''print(1)\nprint(1)'''))
        # nb['cells']=[
        #     nbformat.v4.new_markdown_cell(get_forms(response)[0].getall()),
        #     nbformat.v4.new_markdown_cell(response.css("h1")[0].getall()),
        #     nbformat.v4.new_code_cell('''print(1)\nprint(1)'''),
        #     nbformat.v4.new_code_cell("print(1)")]
        nbformat.write(nb, filename)
        self.log("Saved file {}".format(filename))
