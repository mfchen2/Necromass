from pathlib import Path
from pypdf import PdfReader, PdfWriter
from pypdf.generic import ContentStream, TextStringObject, FloatObject
import pypdfium2 as pdfium

root=Path('/Users/mingfeichen/Manuscript')
source=root/'outputs/manual-20260824-a1/presentations/proteome_metabolome_concordance/assets/proteome_metabolome_concordance.pdf'
writer=PdfWriter(clone_from=source)
page=writer.pages[0]
stream=ContentStream(page.get_contents(),writer)
count=0
for operands,op in stream.operations:
    if op==b'Tm': matrix=operands
    if op==b'Tf': font_name,font_size=operands
    if op not in (b'Tj',b'TJ'): continue
    parts=operands if op==b'Tj' else operands[0]
    text=''.join(str(x) for x in parts if isinstance(x,TextStringObject))
    if text not in ('K','P') and not text.endswith((' / K',' / P')): continue
    old=text[-1]; suffix={'K':'an','P':'h'}[old]
    assert isinstance(parts[-1],TextStringObject) and str(parts[-1]).endswith(old)
    parts[-1]=TextStringObject(str(parts[-1])+suffix)
    font=page['/Resources']['/Font'][font_name].get_object()
    first=int(font['/FirstChar']); widths=font['/Widths']
    extra=sum(float(widths[ord(c)-first]) for c in suffix)/1000*float(font_size)*float(matrix[0])
    if op==b'TJ':
        parts.insert(0,FloatObject(extra/(float(font_size)*float(matrix[0]))*1000))
    else:
        matrix[4]=FloatObject(float(matrix[4])-extra*.5)
    count+=1
assert count==26,count
page.replace_contents(stream)
out=root/'Supplementary_FigureS4_Kan_Ph.pdf'
writer.write(out)
doc=pdfium.PdfDocument(str(out))
doc[0].render(scale=400/72).to_pil().save(root/'Supplementary_FigureS4_Kan_Ph.png')
doc[0].render(scale=1.4).to_pil().save('/tmp/S4_preview.png')
print('Updated labels:',count)
