
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ttest_ind, f_oneway as anova, mannwhitneyu, median_test
from scipy.stats import sem
import pandas as pd
class Annots():
    """The class will generate an errorbar of each group based on the sem
    it will also annotate statistical differences based on a pair group analysis
    It uses only categorical data differences

    plots supported = boxplot, stripplot, violinplot, barplot

    x and hue most be categorical or string"""

    def __init__(self, ax, x, y, hue=None, comp=[], data=None):
        self.x = x
        self.y = y
        self.hue = hue
        self.ax =ax
        self.data = data.astype({self.x:"string", self.hue:"string"})
        self.handles, self.labels = ax.get_legend_handles_labels()
        self.x_labels = ax.get_xticklabels()
        self.x_label_0 = [a.get_text() for a in self.x_labels]
        self.template_labels_0 = self.x_label_0 + self.labels
        self.n_p = int(len(self.labels))
        self.n_x = int(len(self.x_labels))
        self.y=y
        if self.n_p % 2 == 0:
            self.sub_x = 0.8/(2*self.n_p)
            self.n_p_0 = self.n_p - 1
            self.sub_x_p_0 = range(-((self.n_p_0)), (self.n_p_0)+1, 2)
            self.sub_x_p = [i for i in self.sub_x_p_0 if i != 0]
        else:
            self.sub_x = 0.8/self.n_p
            self.sub_x_p = list(range(-(self.n_p-1)//2, (self.n_p//2)+1, 1))
        self.x_position = {self.x_labels[i].get_text():i for i in range(self.n_x)} # xticks coordinates of each label order as in figure
        self.abs_position = {a: {self.labels[i]:self.sub_x_p[i] for i in range(self.n_p)} for a in self.x_label_0} # position of each sub x category as in figure

        self.table_pos = self.data.sort_values([self.x, self.hue], key= lambda x: [self.template_labels_0.index(a) for a in x]).groupby([self.x, self.hue], sort=False)
        self.x_pos_labels = [a[0] for a in self.table_pos.indices.keys()]
        self.x_pos_legends = [a[1] for a in self.table_pos.indices.keys()]
        self.x_pos_0 = []

        self.template_labels = [] # all correctly order possible pairs in the figure  # check if needed
        for i in range(len(self.x_label_0)):
            for j in range(len(self.labels)):
                self.template_labels.append((self.x_label_0[i], self.labels[j]))
        tup2_0 = [str(a) for a in self.x_pos_legends]  #
        tup2 = list(zip(self.x_pos_labels, tup2_0))  # possible pairs as data permitted not order as in the figure
        self.data_or_0 = []
        for a in self.template_labels:
            if a in tup2:
                self.data_or_0.append(a)
        self.data_or = [(self.data_or_0[i], self.data_or_0[i + 1]) for i in range(len(self.data_or_0) - 1)]  # This the possible data groups, ordered correctly
        if comp == []:    # general analysis of each adjacent pair of data groups
            self.comp_empty = True
            self.comp = self.data_or
            #self.comp_fig = [(a[i], a[i+1]) for a in self.all_adjacent for i in range(len(a)-1)] # all possible values in the figure in the figure order
            #self.comp_sorted = sorted(self.comp, key=self.comp_fig)
        else:
            self.comp_empty = False
            tup2_1 = comp
            self.comp_0 = []
            for a in tup2_1:
                for b in self.template_labels:
                    if b in a:
                        self.comp_0.append(b)
            self.comp = [(self.comp_0[i], self.comp_0[i+1]) for i in range(0, len(self.comp_0)-1, 2)]

        self.n_comp = len(self.comp)

        getx1, getx2 = self.ax.get_xlim()
        ax.set_xlim(left=getx1, right=getx2)
        self.yticks = ax.get_yticks()
        self.y_distance = self.yticks[0] - self.yticks[1]

    def max_val(self):
        self.Maxi = []
        self.length_line = []
        self.groups_between = self.data_or
        if self.comp_empty == True:
            Maxi_0 = self.data[
                                (self.data[self.x].astype("string") == self.data_or_0[0][0]) &
                                (self.data[self.hue].astype("string") == self.data_or_0[0][1])
                                ][self.y].max()
            for i in range(1, len(self.data_or_0)):
                self.groups = self.data[
                                        (self.data[self.x].astype("string") == self.data_or_0[i][0]) &
                                        (self.data[self.hue].astype("string") == self.data_or_0[i][1])
                                        ][self.y].max()
                if Maxi_0 < self.groups:
                    self.Maxi.append(self.groups)
                else:
                    self.Maxi.append(Maxi_0)
                Maxi_0 = self.groups
        else:  # Find the groups between the chosen groups for statistical analysis
            for i in range(self.n_comp): #self.n_comp
                self.groups_between = self.template_labels.copy()   # it will contain groups between the two pairs chosen for annotation
                for j in range(2):
                    n_groups_b = len(self.groups_between)
                    for k in range(n_groups_b):
                        if self.groups_between[k] == self.comp[i][j]:
                            if j == 0:
                                del self.groups_between[0:k]
                                break
                            else:
                                del self.groups_between[k+1:]
                                self.Maxi_0 = self.data[
                                                        (self.data[self.x].astype("string") == self.groups_between[0][0]) &
                                                        (self.data[self.hue].astype("string") == self.groups_between[0][1])
                                                        ][self.y].max()
                                for i in range(1, len(self.groups_between)):
                                    self.groups = self.data[
                                                            (self.data[self.x].astype("string") == self.groups_between[i][0]) &
                                                            (self.data[self.hue].astype("string") == self.groups_between[i][1])
                                                            ][self.y].max()
                                    if self.groups > self.Maxi_0:
                                        self.Maxi_0 = self.groups
                                self.Maxi.append(self.Maxi_0)
                                break
        return self.Maxi


    def sorted_length(self):
        pass

    def testing_var(self, test=ttest_ind):
        self.p_values_t = []
        self.groups_stats = []
        for i in range(self.n_comp):
            self.g1 = self.data[(self.data[self.x].astype("string") == self.comp[i][0][0]) & (self.data[self.hue].astype("string") == self.comp[i][0][1])][self.y]
            self.g2 = self.data[(self.data[self.x].astype("string") == self.comp[i][1][0]) & (self.data[self.hue].astype("string") == self.comp[i][1][1])][self.y]
            try:
                self.p_values_t.append(test(self.g1, self.g2)[1])
            except Exception as e:
                print("\033[1;32m One of the groups is empty for statistical testing is empty\n")

        return self.p_values_t

    def p_values(self, test=ttest_ind):
        tests = self.testing_var(test=test)
        self.ast = ""
        self.p_vals = []
        for i in range(self.n_comp):
            self.p_val = tests[i]
            if self.p_val <= 0.05:
                self.ast = "*"
                if self.p_val < 0.01:
                    self.ast += "*"
                    if self.p_val < 0.001:
                        self.ast += "*"
                        if self.p_val < 0.0001:
                            self.ast += "*"
                            self.p_vals.append(self.ast)
                        else:
                            self.p_vals.append(self.ast)
                    else:
                        self.p_vals.append(self.ast)
                else:
                    self.p_vals.append(self.ast)
            else:
                self.p_vals.append("ns")
        return self.p_vals

    def anot_adjacent(self, offset=0.007, yy=None, marks = 10, y_height_distance = 3, significance=False, ext=0, test=ttest_ind):
        """ It annotates on seaborn stripplot, boxplot, barplot"""
        """ y: is the y position for the annotation
         x: is the amount of observations on the X axis
         sub_x is the distance between subcategories in the X axis
         sub_x_p are the relative positions on the Grid 
         
         If you want to check only significant results, you can use significance as True
         
         ext helps when overlapping of annotations"""

        self.ext = abs(self.y_distance)/y_height_distance
        valores_de_y = self.max_val()
        valores_de_significancia = self.p_values(test=test)
        for i in range(self.n_comp):
            if significance == False:
                text = valores_de_significancia[i]
            else:
                if valores_de_significancia[i] == "ns":
                    continue
                else:
                    text = valores_de_significancia[i]
            if yy == None:
                self.yy = valores_de_y[i] + abs(self.y_distance) / y_height_distance + i*self.ext*ext
            else:
                self.yy = yy
            self.y2 = self.yy + abs(self.y_distance) / marks
            #for i in range(len(self.comp)):
            idx1 = str(self.comp[i][0][0]) #x_position[comp[i][0]]
            idx2 = str(self.comp[i][1][0])
            idsubx1 = self.abs_position[idx1][(self.comp[i][0][1])]
            idsubx2 = self.abs_position[idx1][(self.comp[i][1][1])]
            self.x1 = (self.x_position[idx1] + idsubx1 * self.sub_x)  + offset# first element of comparison x_position[idx1]
            self.x2 = (self.x_position[idx2] + idsubx2 * self.sub_x ) - offset# second element of comparison x_position[idx1]
            try:
                plt.plot([self.x1, self.x1, self.x2, self.x2], [self.yy,self.y2, self.y2, self.yy], color="k")
                plt.text((self.x1  + self.x2) * 0.5, self.y2, text, ha="center",va="bottom")
            except:
                continue

    def errorbars(self, trend="_"):
        # fig = plt.gcf()
        # size = fig.get_size_inches() * fig.dpi
        self.means = self.data.groupby([self.x, self.hue], sort=False)[self.y].mean().sort_index(key=lambda x : [self.template_labels_0.index(a) for a in x])
        self.ss = self.data.groupby([self.x, self.hue], sort=False)[self.y].sem().sort_index(key=lambda x : [self.template_labels_0.index(a) for a in x])
        self.x_pos_0 = [a[0] for a in self.means.index]
        self.x_pos_1 = [a[1] for a in self.means.index]
        self.x_pos_2 = []
        for i in range(len(self.x_pos_0)):
            a = self.x_pos_0[i]
            b = str(self.x_pos_1[i])
            exes = self.x_position[a] + self.abs_position[a][b] * self.sub_x
            self.x_pos_2.append(exes)
        plt.errorbar(self.x_pos_2, self.means, self.ss,  markersize=15, capsize=5, color="k", zorder=1, fmt=trend)


