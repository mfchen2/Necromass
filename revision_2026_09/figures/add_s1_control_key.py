"""Compose the complete original S1 with aligned vector legend keys."""
from pathlib import Path
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
import pypdfium2 as pdfium

root=Path('/Users/mingfeichen/Manuscript')
source=Path('/Users/mingfeichen/Documents/Laptop_Necromass/Necromass_Manuscript/Submission/Necromass_Supplementary_FigureS1.jpg')
image=ImageReader(str(source)); iw,ih=image.getSize()
width=656; height=width*ih/iw
out=root/'Supplementary_FigureS1_control_key.pdf'
c=canvas.Canvas(str(out),pagesize=(710,height+12))
c.drawImage(image,5,6,width=width,height=height)
c.setFont('Helvetica',10)
for fraction in (735/1893,1650/1893):
    y=6+height*(1-fraction)
    c.setStrokeColorRGB(0,0,0); c.setLineWidth(1.5); c.setDash(4,2.5)
    c.line(567,y,582,y)
    c.setDash()
    c.drawString(591,y-3.4,'Untreated control')
c.save()
d=pdfium.PdfDocument(str(out))
d[0].render(scale=400/72).to_pil().save(root/'Supplementary_FigureS1_control_key.png')
d[0].render(scale=1).to_pil().save('/tmp/s1_key_preview.png')
