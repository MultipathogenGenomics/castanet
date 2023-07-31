import time
import pandas as pd
from datetime import datetime
from functools import partial

import reportlab.lib.colors as clr
from reportlab.platypus import Paragraph, Image, Frame, BaseDocTemplate, PageTemplate, HRFlowable, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.lib import pagesizes

from app.utils.pdf_mappings import get_bullet_style, get_fluff_style, get_header_style


class GenerateReport:
    def __init__(self, contigs_fig_path, consensus_fig_path, payload, t_stats, c_stats, m_stats) -> None:
        self.a = payload
        self.t_stats = t_stats   # Target consensuses
        self.c_stats = c_stats   # Remapped consensus
        self.m_stats = m_stats  # Mash between consensus types
        self.depth_df = pd.read_csv(
            f"{self.a['folder_stem']}{self.a['SeqName']}_depth.csv")
        self.depth = self.depth_df[self.depth_df["probetype"]
                                   == self.a['GtOrg']]
        self.contigs_fig_path = contigs_fig_path
        self.consensus_fig_path = consensus_fig_path
        self.depth_fig_path = f"{self.a['folder_stem']}Depth_output/{self.a['GtOrg']}-{self.a['SeqName']}.png"
        self.page_size = pagesizes.portrait(pagesizes.A4)
        self.now = datetime.now()
        self.styles = getSampleStyleSheet()
        self.fname = f"{self.a['folder_stem']}evaluation/{self.a['SeqName']}_{self.a['GtOrg']}_run_report.pdf"

    def footer(self, canvas, doc, content) -> None:
        '''Render content on footer'''
        canvas.saveState()
        w, h = content.wrap(doc.width, doc.bottomMargin)
        content.drawOn(canvas, doc.leftMargin, h)
        canvas.restoreState()

    def header(self, canvas, doc, content) -> None:
        '''Render content on header'''
        canvas.saveState()
        w, h = content.wrap(doc.width, doc.bottomMargin)
        content.drawOn(canvas, doc.leftMargin,
                       doc.height + 0.3 * doc.topMargin)
        content2 = Paragraph(
            f"Report Date: {datetime.now().strftime('%d-%m-%Y')}", get_header_style(self.styles))
        w, h = content2.wrap(doc.width, doc.bottomMargin)
        content2.drawOn(canvas, doc.leftMargin,
                        doc.height + 0.6 * doc.topMargin)
        content3 = Paragraph(
            f"Report Time: {datetime.now().isoformat(timespec='seconds')[-8:]} UTC", get_header_style(self.styles))
        w, h = content3.wrap(doc.width, doc.bottomMargin)
        content3.drawOn(canvas, doc.leftMargin,
                        doc.height + 0.5 * doc.topMargin)
        canvas.restoreState()

    def header_and_footer(self, canvas, doc, header_content, footer_content) -> None:
        '''Call header, footer and page divider render while inside partial function'''
        self.footer(canvas, doc, header_content)
        self.header(canvas, doc, footer_content)
        canvas.line(doc.leftMargin, 100, doc.width, 100)

    def build_logo(self, fpath, l=540, w=90):
        '''Build image object'''
        logo = Image(fpath, l, w)
        logo.hAlign = "CENTER"
        return logo

    def split_stats(self):
        table_data = [self.m_stats.columns[:,].values.astype(
            str).tolist()] + self.m_stats.values.tolist()
        table_data = [i[1:5] for i in table_data]
        vs_genome_data, vs_gs_data = [], []
        for i in range(len(table_data)):
            if not i == 0:
                table_data[i][0] = table_data[i][0].split("/")[-1]
            if i == 0:
                vs_genome_data.append(
                    table_data[i]), vs_gs_data.append(table_data[i])
            elif i <= 3:
                vs_genome_data.append(table_data[i])
            elif i > 3:
                vs_gs_data.append(table_data[i])
        return vs_genome_data, vs_gs_data

    def get_summary_stats(self):
        stats = [["run_time", "n_ref_org_reads", "mean_depth", "mean_depth_std"]]
        stats.append([
            f"{round((time.time() - self.a['StartTime']) / 60)} min",
            f"{int(self.depth['n_reads_all'].iloc[0])}",
            f"{round(self.depth['depth_mean'].iloc[0])}",
            f"{round(self.depth['depth_std'].iloc[0])}"
        ])
        return stats

    def build_table(self, data, dims=(10*cm, 2.5*cm, 1.0*cm, 3.5*cm)):
        '''Construct table style and dimensions, insert data'''
        style = []
        for i, _ in enumerate(data):
            if i % 2 != 0:
                style.append(
                    ('BACKGROUND', (0, i), (-1, i), clr.Color(0, 183, 155, 0.1)))
        style = TableStyle(style)
        tbl = Table(data, colWidths=dims)
        tbl.setStyle(style)
        return tbl, style

    def main(self) -> str:
        '''Entrypoint: build canvas and call renderers'''
        doc = BaseDocTemplate(self.fname,
                              pagesize=self.page_size,
                              leftMargin=1 * cm,
                              rightMargin=1 * cm,
                              topMargin=5 * cm,
                              bottomMargin=0.1 * cm)
        bullet_point_style = get_bullet_style(self.styles)
        story = []

        '''Header/footer static content'''
        footer_content = Paragraph(f"Sample: {self.a['ExpName']}, Ref Organism: {self.a['GtOrg']}", get_fluff_style(
            self.styles))  # Image("./app/assets/osprey.png", 122/2, 142/2)
        header_content = Paragraph(
            "", get_fluff_style(self.styles))

        padding = dict(
            leftPadding=35,
            rightPadding=35,
            topPadding=120,
            bottomPadding=120)

        '''Create a Frame for the Flowables (Paragraphs and such)'''
        frame = Frame(0, 0, *A4, **padding)

        '''Add the Frame to the template and tell the template to call draw static for each page'''
        template = PageTemplate(id='test', frames=[frame], onPage=partial(
            self.header_and_footer, header_content=header_content, footer_content=footer_content))

        '''Add the template to the doc'''
        doc.addPageTemplates([template])

        '''Separate stats for tables'''
        vs_genome_data, vs_gs_data = self.split_stats()
        summary_data = self.get_summary_stats()

        '''Title'''
        title_style = self.styles['Heading1']
        title_style.alignment = 1
        title_style.fontSize = 30
        whitesp = Paragraph("&nbsp;", self.styles["Normal"])
        title = Paragraph(
            f"Castanet Eval Report", title_style)
        story.append(title)
        story.append(whitesp)

        '''Header HLine'''
        hline = HRFlowable(width="100%", thickness=1, lineCap='round', color=clr.black,
                           spaceBefore=1, spaceAfter=1, hAlign='CENTER', vAlign='BOTTOM', dash=[10, 10])
        story.append(hline)

        '''Run stats details'''
        story.append(Paragraph("Run Summary", self.styles["Heading2"]))
        '''Bullet points for run stats'''
        summary_table, summary_tbl_style = self.build_table(
            summary_data, dims=(3.0*cm, 3.0*cm, 3.0*cm, 3.0*cm))
        summary_table.setStyle(summary_tbl_style)
        story.append(summary_table)
        '''Depth Image'''
        story.append(self.build_logo(self.depth_fig_path, 240, 180))
        story.append(Paragraph("&nbsp;", self.styles["Normal"]))
        story.append(Paragraph("&nbsp;", self.styles["Normal"]))
        story.append(Paragraph("&nbsp;", self.styles["Normal"]))

        '''Remapped consensus Key Stats'''
        story.append(Paragraph(
            "Castanet (Re-mapped) Consensus Statistics", self.styles["Heading2"]))
        '''Vs Genome table'''
        vs_gs_table, genome_tbl_style = self.build_table(
            self.c_stats, (3.0*cm, 3.0*cm, 3.0*cm, 3.0*cm, 3.0*cm))
        vs_gs_table.setStyle(genome_tbl_style)
        story.append(vs_gs_table)
        story.append(Paragraph("&nbsp;", self.styles["Normal"]))

        '''Consensus vs known Details'''
        story.append(
            Paragraph("All consensus vs Known Viral Sequence", self.styles["Heading2"]))
        '''Vs Genome table'''
        genome_tbl, genome_tbl_style = self.build_table(vs_genome_data)
        genome_tbl.setStyle(genome_tbl_style)
        story.append(genome_tbl)
        story.append(Paragraph("&nbsp;", self.styles["Normal"]))
        story.append(Paragraph("&nbsp;", self.styles["Normal"]))
        '''Consensus Image'''
        story.append(Paragraph(
            "Graph shows similarity of consensuses and known viral GenBank sequence", self.styles["Normal"]))
        story.append(self.build_logo(self.consensus_fig_path))
        story.append(Paragraph("&nbsp;", self.styles["Normal"]))
        story.append(Paragraph("&nbsp;", self.styles["Normal"]))
        story.append(Paragraph("&nbsp;", self.styles["Normal"]))

        '''Consensus vs gold standard details'''
        story.append(Paragraph(
            "Castanet consensus vs Gold Standard Consensus", self.styles["Heading2"]))
        '''Vs Genome table'''
        vs_gs_table, genome_tbl_style = self.build_table(vs_gs_data)
        vs_gs_table.setStyle(genome_tbl_style)
        story.append(vs_gs_table)
        story.append(Paragraph("&nbsp;", self.styles["Normal"]))

        '''Contigs details'''
        story.append(Paragraph(
            "Contigs (individual target consensuses) vs Targets", self.styles["Heading2"]))
        '''Contig coverage table'''
        genome_tbl, genome_tbl_style = self.build_table(self.t_stats, dims=(
            6.25*cm, 1.75*cm, 1.75*cm, 1.75*cm, 1.75*cm, 1.75*cm, 1.75*cm, 1.75*cm, 1.75*cm))
        genome_tbl.setStyle(genome_tbl_style)
        story.append(genome_tbl)
        story.append(Paragraph("&nbsp;", self.styles["Normal"]))
        story.append(Paragraph("&nbsp;", self.styles["Normal"]))
        story.append(Paragraph(
            "Graph shows similarity of all relevant targets, with target consensuses stacked on top", self.styles["Normal"]))
        '''Contigs Image'''
        story.append(self.build_logo(self.contigs_fig_path))
        story.append(Paragraph("&nbsp;", self.styles["Normal"]))
        story.append(Paragraph("&nbsp;", self.styles["Normal"]))

        doc.build(story)
        return self.fname
