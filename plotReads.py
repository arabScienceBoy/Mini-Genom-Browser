from __future__ import division
import re
import Tkinter as t
from tkFileDialog import askopenfilename
from tkMessageBox import showerror
from Bio import SeqIO, AlignIO, SearchIO


import pysam
import canvasvg as svg

# todo: mach das Limit der reads pro Spalte auch wenn die Reads teilweise oder
# vollstaendig nicht mehr sichtbar!

class DrawGenome(object):

    def __init__(self, consStart, consEnd, genomDNA, reads, posScale=1,
                 genomScale=1, width=900, height=700):
        self.consStart = consStart
        self.consEnd = consEnd
        self.consBlockSize = 50  # 500
        self.blockSize = 500  # self.consBlockSize
        self.genomStart = self.consStart
        self.genomEnd = self.blockSize
        self.genomDNA = genomDNA
        self.posScale = posScale
        self.genomScale = genomScale
        self.width = width
        self.height = height
        self.reads = reads

        root = t.Tk()
        root.title("Mapping Reads to Reference Genome")
        self.initGui(root)
        root.eval('tk::PlaceWindow %s center' % root.winfo_pathname(root.winfo_id()))
        t.mainloop()

    def reset(self):
        self.genomStart = self.consStart
        self.genomEnd = self.consEnd

    def drawAll(self):
        self.drawGenom()
        self.drawReads()

    def initGui(self, root):
        frame = t.Frame(root, borderwidth=6, relief='raised', width=self.width)
        frame.pack(fill='both', padx=10)

        l = t.Label(frame, text="control panel")

        l.grid(row=0)
        l.pack(fill='both')

        button1 = t.Button(frame, text="zoom out", command=self.zoomOut)
        button1.pack(side="top")

        button2 = t.Button(frame, text="zoom in", command=self.zoomIn)
        button2.pack(side="top")

        button3 = t.Button(frame, text="move forward", command=self.moveRight)
        button3.pack(side="right")

        button4 = t.Button(frame, text="move backward", command=self.moveLeft)
        button4.pack(side="left")

        positionEntry1 = t.Entry(frame, width=15)
        positionEntry2 = t.Entry(frame, width=15)
        positionEntry1.insert(0, 'hatem')
        positionEntry2.insert(9, 'Shugaa')
        positionEntry1.pack(side="left")
        positionEntry2.pack(side="right")

        self.convas = t.Canvas(root, width=self.width, height=self.height,
                               borderwidth=5, relief='raised')

        self.convas.pack(fill='both', padx=10)

        self.drawGenom()

    def drawGenom(self):
        # problem: wenn man zoom out, mit consBlockSize = 500, Ende:40000
        # wird verkleinert bis 450000 und nicht bis 40000!
        if self.genomStart + self.blockSize >= self.consEnd:
            last_blockSize = self.blockSize
            now_blockSize = abs(self.consEnd - self.genomStart)
            self.genomStart -= abs(now_blockSize - last_blockSize)

        if self.genomStart < self.consStart:
            self.genomStart = self.consStart

        genomLineEnd = self.width - 15
        genomLineStart = 30

        self.convas.create_line(genomLineStart, 50,
                                genomLineEnd, 50, width=20, fill='green')

        if len(self.reads) > 1:
            self.convas.create_line(genomLineStart, 270,
                                    genomLineEnd, 270, width=2, fill='black')
            self.convas.create_text(genomLineStart - 13 + 3.5, 150, text='B1')
            self.convas.create_text(genomLineStart - 13 + 3.5, 380, text='B2')

        # draw genome axis
        lineLen = self.posScale

        toAddScaled = self.blockSize / self.posScale
        toAddReal = (genomLineEnd - 30) / self.posScale

        pos_scaled = self.genomStart
        pos_real = genomLineStart - 4

        self.convas.create_text(pos_real + 3.5, 20, text='%d' % pos_scaled)
        self.convas.create_line(pos_real + 3.5, 500, pos_real + 3.5, 30)

        for _ in range(lineLen):
            pos_scaled += toAddScaled
            pos_real += toAddReal
            self.convas.create_text(pos_real + 3.5, 20, text='%d' % pos_scaled)
            self.convas.create_line(pos_real + 3.5, 500, pos_real + 3.5, 30)

        self.convas.pack()
        self.drawReads()

    def drawReads(self):
        genomLineEnd = self.width - 15
        rSep = 0  # 30 Trennung der Reads
        after = 0
        count = 0

        returnToTop = 0
        numReads = len(self.reads)
        for r in xrange(0, numReads):
            # recHigRead ist der Read, wo noch damit viele ueberlappungen gibt!
            # recHigRead = 0: Read mit index 0 hat noch ueberlappungen!
            recHigRead = 0
            for i in xrange(0, len(self.reads[r])):
                rStart = self.reads[r][i].get_blocks()[0][0]
                rEnd = self.reads[r][i].get_blocks()[0][1]
                rrEnd = self.reads[r][recHigRead].get_blocks()[0][1]

                matches = self.isMismatch(self.reads[r][i])
                color = ''
                if matches == 0:
                    color = 'green'
                elif matches == 1:
                    color = 'orange'
                elif matches == 2:
                    color = 'orange red'
                else:
                    color = 'red'

                # if the reads are outside the display area from the right, don't
                # draw
                if rEnd > self.genomEnd:
                    # if the whole read is outside the display area, don't show
                    # anything!
                    if rStart > self.genomEnd:
                        if (rrEnd + rSep) < rStart:
                            recHigRead = i
                            after = returnToTop
                            break
                        else:
                            after += 6.5
                            count += 1
                        # continue

                    # if only part of the read is showed, then cut it and show the part
                    # that can be showed!
                    else:
                        readStart = ((rStart - self.genomStart) /
                                     self.blockSize) * (genomLineEnd - 30)
                        readEnd = ((rEnd - self.genomStart -
                                    abs(self.genomEnd - rEnd)) / self.blockSize) * (genomLineEnd - 30)
                        # test if there is overlapp with the recent read
                        if (rrEnd + rSep) < rStart:
                            recHigRead = i
                            after = returnToTop
                            count = 0
                            self.convas.create_line(readStart + 31,
                                                    after + 63, readEnd + 30,
                                                    after + 63, width=5, fill=color)
                            after += 6.5
                        else:
                            if count > 30:
                                continue
                            else:
                                self.convas.create_line(readStart + 31,
                                                        after + 63, readEnd + 30,
                                                        after + 63, width=5, fill=color)
                                after += 6.5
                                count += 1
                                # continue

                # if the reads are outside the display area from the left, don't draw
                # the outside reads! fix the problem with svg!
                elif rStart < self.genomStart:
                    # if the whole read is outside the display area, don't show
                    # anything!
                    if rEnd < self.genomStart:
                        if (rrEnd + rSep) < rStart:
                            recHigRead = i
                            after = returnToTop
                            after += 6.5
                            count = 0
                        else:
                            after += 6.5
                            count += 1
                        # continue

                    # if only part of the read is showed, then cut it and show the
                    # part that can be showed!
                    else:
                        readStart = ((rStart - self.genomStart +
                                      abs(rStart - self.genomStart)) /
                                     self.blockSize) * (genomLineEnd - 30)
                        readEnd = ((rEnd - self.genomStart) /
                                   self.blockSize) * (genomLineEnd - 30)

                        if (rrEnd + rSep) < rStart:
                            recHigRead = i
                            after = returnToTop
                            count = 0
                            self.convas.create_line(readStart + 31,
                                                    after + 63, readEnd + 30,
                                                    after + 63, width=5, fill=color)
                            after += 6.5

                        else:
                            if count > 30:
                                continue
                            else:
                                self.convas.create_line(readStart + 31,
                                                        after + 63, readEnd + 30,
                                                        after + 63, width=5, fill=color)
                                after += 6.5
                                count += 1
                                # continue

                # problem: see 12100
                # the problem is because of the position of "after"!
                # should "after" come after recHigRead?
                else:
                    readStart = ((rStart - self.genomStart) /
                                 self.blockSize) * (genomLineEnd - 30)
                    readEnd = ((rEnd - self.genomStart) /
                               self.blockSize) * (genomLineEnd - 30)

                    if (rrEnd + rSep) < rStart:
                        recHigRead = i
                        after = returnToTop
                        count = 0
                        self.convas.create_line(readStart + 31,
                                                after + 63, readEnd + 30,
                                                after + 63, width=5, fill=color)
                        after += 6.5

                    else:
                        # show only 30 reads pro Spalte (wichtig wegen zwei
                        # Bams)
                        if count > 30:
                            continue
                        else:
                            self.convas.create_line(readStart + 31,
                                                    after + 63, readEnd + 30,
                                                    after + 63, width=5, fill=color)
                            after += 6.5
                            count += 1
            returnToTop = 210
            count = 0
            after = returnToTop

        self.convas.pack()
        svg.saveall('probe.svg', self.convas)

    def zoomOut(self):
        if self.consEnd <= self.genomEnd:
            if self.consStart >= self.genomStart:
                print "genome is maximal zoomed out!"
                return
            else:
                self.genomStart -= self.consBlockSize
                self.blockSize += self.consBlockSize
        elif self.consStart >= self.genomStart:
            if self.consEnd <= self.genomEnd:
                print "genome is maximal zoomed out!"
                return
            else:
                self.blockSize += self.consBlockSize
                self.genomEnd += self.consBlockSize
        else:
            self.genomStart -= self.consBlockSize
            self.genomEnd += self.consBlockSize
            self.blockSize += 2 * self.consBlockSize

        self.convas.delete('all')
        self.drawGenom()
        # self.drawReads()

    def zoomIn(self):
        # Problem: wenn man weiter vergroessert, kommt sowas:
        # Start: 341, End: 334
        # End kleiner als Start! das darf nie passieren!
        if self.blockSize <= self.consBlockSize:
            self.blockSize = self.consBlockSize
            return
        if self.consEnd <= self.genomEnd:
            if self.consStart >= self.genomStart:
                print "genome is maximal zoomed out!"
                self.genomStart += self.consBlockSize
                self.genomEnd -= self.consBlockSize
                self.blockSize -= self.consBlockSize
                return
            else:
                self.genomStart += self.consBlockSize
                self.blockSize -= self.consBlockSize
        elif self.consStart >= self.genomStart:
            if self.consEnd <= self.genomEnd:
                print "genome is maximal zoomed in!"
                #self.genomEnd = self.consEnd
                return
            else:
                self.genomEnd -= self.consBlockSize
                self.blockSize -= self.consBlockSize
        else:
            if self.blockSize == 2 * self.consBlockSize:
                self.genomEnd -= self.consBlockSize
                self.blockSize = self.consBlockSize
            else:
                self.genomStart += self.consBlockSize
                self.genomEnd -= self.consBlockSize
                self.blockSize -= 2 * self.consBlockSize

        self.convas.delete('all')
        self.drawGenom()

    def moveRight(self):
        if self.consEnd <= self.genomEnd:
            print "end of genome........!"
            self.genomEnd = self.consEnd
            return
        else:
            if self.blockSize == 2 * self.consBlockSize:
                self.genomEnd += self.consBlockSize  # self.blockSize
            else:
                self.genomStart += self.consBlockSize  # self.blockSize
                self.genomEnd += self.consBlockSize  # self.blockSize

        self.convas.delete('all')
        self.drawGenom()

    def moveLeft(self):
        if self.consStart >= self.genomStart:
            print "end of genome........!"
            self.genomStart = self.consStart
            return
        else:
            self.genomStart -= self.consBlockSize  # self.blockSize
            self.genomEnd -= self.consBlockSize  # self.blockSize

        self.convas.delete('all')
        self.drawGenom()

    def scaleGenom(self, tt):
        self.genomEnd = ((int(tt) / 100) * self.consEnd) + self.genomStart
        if int(tt) == 0:
            self.blockSize = 500
        else:
            self.blockSize = (int(tt) / 100) * self.consEnd

        if self.genomStart > 0:
            self.scal2.config(from_=self.genomEnd)
            # self.scal2.set(self.genomEnd)
            # pass
        else:
            pass
            # self.scal2.set(self.genomEnd)
            # self.scal2.config(from_=self.genomEnd)
            #

        self.convas.delete('all')
        self.drawGenom()

    def moveInGenom(self, tt):
        self.genomEnd = int(tt)
        self.genomStart = abs(self.genomEnd - self.blockSize)

        self.convas.delete('all')
        self.drawGenom()

    def isMismatch(self, read):
        reg = '([0-9]+|[A-Z][0-9]+)'
        md = (read.tags[1][1])
        """if re.match('[0-9]+\s*$', md):
            return 0
        else:
            return 1 """

        result = re.findall(reg, md)
        if len(result) == 1:
            if str.isalpha(result[0][0]):
                return 1
            else:
                return 0
        elif len(result) == 2:
            return 1
        elif len(result) == 3:
            return 2
        else:
            return 3


genBankPath = ''
bamFile1 = ''
bamFile2 = ''


def loadFiles():
    window = t.Tk()
    window.title("Load GenomDB and Bam Files")
    window.rowconfigure(5, weight=1)
    window.columnconfigure(5, weight=1)

    buttonGeDB = t.Button(text="choose GenomDB",
                          command=load_GenBank,  width=20)
    buttonGeDB.grid(row=1, column=1)

    buttonBam1 = t.Button(text="choose Bam file 1",
                          command=load_Bam1, width=20)
    buttonBam1.grid(row=2, column=0)

    buttonBam2 = t.Button(text="choose Bam file 2",
                          command=load_Bam2, width=20)
    buttonBam2.grid(row=2, column=2)

    buttonOk = t.Button(text="Press to Continue",
                        command=window.destroy, width=20)
    buttonOk.grid(row=3, column=1)
    
    window.eval('tk::PlaceWindow %s center' % window.winfo_pathname(window.winfo_id()))
    t.mainloop()


def load_GenBank():
    fileName = askopenfilename(filetypes=(("GenBank", "*.gbk"),
                                          ("All files", "*.*")))
    global genBankPath
    if fileName:
        try:
            genBankPath = fileName
        except:
            showerror("Open Source File",
                      "Failed to read file\n'%s'" % fileName)
        return


def load_Bam1():
    fileName = askopenfilename(filetypes=(("Bam file", "*.bam"),
                                          ("All files", "*.*")))
    global bamFile1
    if fileName:
        try:
            bamFile1 = fileName
        except:
            showerror("Open Source File",
                      "Failed to read file\n'%s'" % fileName)
        return


def load_Bam2():
    fileName = askopenfilename(filetypes=(("Bam file", "*.bam"),
                                          ("All files", "*.*")))
    global bamFile2
    if fileName:
        try:
            bamFile2 = fileName
        except:
            showerror("Open Source File",
                      "Failed to read file\n'%s'" % fileName)
        return


def initGenomBams(genBankPat, bams):
    genomRecords = SeqIO.parse(genBankPat, "genbank")
    genom = next(genomRecords)

    genomStart = 0
    genomEnd = len(genom.seq)
    genomText = genom.seq

    bamfiles = [pysam.AlignmentFile(bamPath, "rb") for bamPath in bams]

    readsList1 = []
    readsList2 = []
    reads = []

    if len(bamfiles) == 2:
        reads_1 = bamfiles[0].fetch()
        for r1 in reads_1:
            readsList1 += [r1]

        reads_2 = bamfiles[1].fetch()
        for r2 in reads_2:
            readsList2 += [r2]

        reads = [readsList1, readsList2]
        DrawGenome(genomStart, genomEnd, genomText, reads, 10, 4)

    else:
        reads_1 = bamfiles[0].fetch()
        for r1 in reads_1:
            readsList1 += [r1]
        reads = [readsList1]
        DrawGenome(genomStart, genomEnd, genomText, reads, 10, 4)


loadFiles()


if genBankPath != '' and bamFile1 != '' and bamFile2 != '':
    initGenomBams(genBankPath, [bamFile1, bamFile2])
elif genBankPath and bamFile1:
    initGenomBams(genBankPath, [bamFile1])
elif genBankPath and bamFile2:
    initGenomBams(genBankPath, [bamFile2])
else:
    showerror("Error:", "You didn't choose genbank file and bam files...!")
