from pathlib import Path
from typing import List
from typing import Any, Dict, Union
from sphinx.config import Config
from sphinx.application import Sphinx
from docutils import nodes
from docutils.nodes import Node
from docutils.parsers.rst import directives
from sphinx.locale import _
from sphinx.util.docutils import SphinxDirective
from sphinx.util.fileutil import copy_asset
from sphinx.util import logging

logger = logging.getLogger(__name__)


class enumerable(nodes.Admonition, nodes.Element):
    pass


def EnumerableFactory(enumerable_type, directivetype):

    class Enumerable(SphinxDirective):

        typ = directivetype
        required_arguments = 0
        optional_arguments = 1
        final_argument_whitespace = True
        has_content = True

        option_spec = {
            'name': directives.unchanged,
            'class': directives.class_option,
        }

        def run(self) -> List[Node]:

            self.assert_has_content()
            env = self.env

            name = self.options.get("name", "")

            if name:
                self.options["noindex"] = False
                node_id = name
            else:
                self.options["noindex"] = True
                node_id = env.new_serialno()

            ids = [node_id]

            classes, class_name = ['enumerable'], self.options.get("class", [])
            if class_name:
                classes.extend(class_name)

            title = ""
            if self.arguments != []:
                title = "  ({})".format(self.arguments[0])

            the_node = enumerable_type(rawsource='\n'.join(self.content),
                                       classes=classes,
                                       **self.options)
            the_node += nodes.title(_(title), _(title))

            the_node["ids"].extend(ids)
            the_node["type"] = self.typ

            self.state.nested_parse(self.content,
                                    self.content_offset, the_node)
            self.add_name(the_node)

            return [the_node]

    return Enumerable


def visit_the_node(self, node: Node) -> None:
    self.visit_admonition(node)


def depart_the_node(self, node: Node) -> None:
    self.depart_admonition(node)


def init_numfig(app: Sphinx, config: Config) -> None:
    """Initialize proof numfig format."""
    config["numfig"] = True
    numfig_format = dict()
    for en in config['enumerable_envs']:
        numfig_format[en[0]] = "{} %s".format(en[1])
    config.numfig_format = numfig_format


def copy_asset_files(app: Sphinx, exc: Union[bool, Exception]) -> None:
    static_path = (
        Path(__file__).parent.parent.joinpath("static",
                                              "enumerable.css").absolute()
    )
    asset_files = [str(static_path)]

    if exc is None:
        for path in asset_files:
            copy_asset(path,
                       str(Path(app.outdir).joinpath("_static").absolute()))


def new_enumerable(
        app: Sphinx, envname: str, displayname: str, counter: str) -> None:

    classname = 'enumerable_{}'.format(envname)
    enumerable_node = type(classname, (enumerable,), {})
    globals()[classname] = enumerable_node  # important for pickling

    app.add_directive(envname, EnumerableFactory(enumerable_node, displayname))
    app.add_enumerable_node(enumerable_node, counter, None,
                            html=(visit_the_node, depart_the_node),
                            latex=(visit_the_node, depart_the_node),
                            text=(visit_the_node, depart_the_node))


def install_extension(app: Sphinx, config: Config) -> None:

    for en in config['enumerable_envs']:
        new_enumerable(app, *en)


def setup(app: Sphinx) -> Dict[str, Any]:

    enum_envs = [
        ('acknowledgment', 'Acknowledgment', 'acknowledgment'),
        ('assertion', 'Assertion', 'assertion'),
        ('assumption', 'Assumption', 'assumption'),
        ('axiom', 'Axiom', 'axiom'),
        ('case', 'Case', 'case'),
        ('claim', 'Claim', 'claim'),
        ('conclusion', 'Conclusion', 'conclusion'),
        ('condition', 'Condition', 'condition'),
        ('conjecture', 'Conjecture', 'conjecture'),
        ('corollary', 'Corollary', 'corollary'),
        ('criterion', 'Criterion', 'criterion'),
        ('definition', 'Definition', 'definition'),
        ('example', 'Example', 'example'),
        ('exercise', 'Exercise', 'exercise'),
        ('hypothesis', 'Hypothesis', 'hypothesis'),
        ('lemma', 'Lemma', 'lemma'),
        ('notation', 'Notation', 'notation'),
        ('observation', 'Observation', 'observation'),
        ('problem', 'Problem', 'problem'),
        ('property', 'Property', 'property'),
        ('proposition', 'Proposition', 'proposition'),
        ('question', 'Question', 'question'),
        ('remark', 'Remark', 'remark'),
        ('summary', 'Summary', 'summary'),
        ('theorem', 'Theorem', 'theorem'),
    ]

    app.add_config_value('enumerable_envs', enum_envs, 'env')

    app.connect("build-finished", copy_asset_files)
    app.connect("config-inited", install_extension)
    app.connect("config-inited", init_numfig)

    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }


class EnumerableException(Exception):
    pass


for name in globals().copy():
    if name.startswith('enumerable_'):
        raise EnumerableException(
            f"Enumerable Internal Error: '{name}' in globals()")
