"""Rebuild supplementary pages with uncropped images and aligned captions."""
import ast
from pathlib import Path
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph
from reportlab.lib.utils import ImageReader
from pypdf import PdfReader

BASE=Path('/Users/mingfeichen/Manuscript')
ROOT=Path('/Users/mingfeichen/Documents/Laptop_Necromass/Necromass_Manuscript/Submission')
OUT=BASE/'Supplementary_Figures_Combined_corrected.pdf'
tree=ast.parse((BASE/'make_supplementary_figures_pdf.py').read_text())
captions=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign)
              and any(isinstance(t,ast.Name) and t.id=='CAPTIONS' for t in n.targets))
captions['S5']+=' The plotted term "adsorption" denotes operationally measured aqueous-phase removal; loss from the supernatant does not distinguish surface adsorption from other removal processes.'
captions['S4']+=' Al, aluminum; Kan, kanamycin; Ph, phage.'
style=ParagraphStyle('caption',fontName='Helvetica',fontSize=10,leading=13)
w,h=A4
margin=31
c=canvas.Canvas(str(OUT),pagesize=A4)
c.setTitle('Supplementary Figures S1-S5')
for i in range(1,6):
    label=f'S{i}'
    path=ROOT/f'Necromass_Supplementary_Figure{label}.{ "jpg" if i==1 else "png"}'
    if i==1: path=BASE/'Supplementary_FigureS1_control_key.png'
    if i==3: path=BASE/'Supplementary_FigureS3_Kan_Ph.png'
    if i==4: path=BASE/'Supplementary_FigureS4_Kan_Ph.png'
    p=Paragraph(f'<b>Supplementary Figure {label}.</b> '+captions[label],style)
    _,ph=p.wrap(w-2*margin,h)
    caption_y=44
    image_bottom=caption_y+ph+14
    image_top=h-55
    image=ImageReader(str(path)); iw,ih=image.getSize()
    scale=min((w-2*margin)/iw,(image_top-image_bottom)/ih)
    dw,dh=iw*scale,ih*scale
    x=(w-dw)/2; y=image_top-dh
    assert y>=image_bottom-1e-6 and y+dh<=h-55+1e-6
    c.setFont('Helvetica-Bold',12)
    c.drawString(margin,h-32,f'Supplementary Figure {label}')
    c.drawImage(image,x,y,width=dw,height=dh,mask='auto')
    p.drawOn(c,margin,caption_y)
    c.setFont('Helvetica',8)
    c.drawCentredString(w/2,19,f'Supplementary Information | {i}')
    c.showPage()
c.save()
r=PdfReader(OUT)
assert len(r.pages)==5
assert 'kanamycin (Kan) or phage (Ph)' in r.pages[2].extract_text()
assert 'does not distinguish surface adsorption' in r.pages[4].extract_text()
print(OUT)
