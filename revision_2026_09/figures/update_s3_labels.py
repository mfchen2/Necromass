from pypdf import PdfReader, PdfWriter
from pypdf.generic import ContentStream, TextStringObject
from pathlib import Path

root=Path('/Users/mingfeichen/Manuscript')
r=PdfReader('/Users/mingfeichen/PCoA_Gower_untargeted_metabolites_012926.pdf')
p=r.pages[0]
stream=ContentStream(p.get_contents(),r)
changed=[]
def replace(value):
    if isinstance(value,TextStringObject) and str(value) in ('K','P'):
        changed.append(str(value))
        return TextStringObject({'K':'Kan','P':'Ph'}[str(value)])
    return value
for operands,operator in stream.operations:
    if operator==b'Tj': operands[0]=replace(operands[0])
    elif operator==b'TJ':
        for i,value in enumerate(operands[0]): operands[0][i]=replace(value)
assert sorted(changed)==['K','P'],changed
p.replace_contents(stream)
w=PdfWriter(); w.add_page(p)
w.write(root/'Supplementary_FigureS3_Kan_Ph.pdf')
print(changed)
