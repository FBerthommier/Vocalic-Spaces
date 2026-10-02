"""Strict verification of the open-access edition (paper/paper.tex) against
the original submission source (Paper.tex of the JASA supplement).

Checks, after normalizing whitespace and comments:

1. word-level diff of the whole body -- every difference must belong to the
   documented intentional changes (template swap: footnotes pointing to this
   repository, plainnat bibliography, caption package, centered figures,
   dropped submission-only end matter);
2. equation environments (equation/align/subequations/empheq) -- must match
   exactly;
3. table rows (cell by cell);
4. figure captions;
5. multimedia (Mm.) descriptions;
6. \\textipa{...} phonetic strings (multiset comparison; the original has
   each caption twice because of the submission-format caption list).

Usage:  python verify_paper.py [path/to/original/Paper.tex]
        (default: ../../matlab_supplement-adjacent original location,
         or any copy of the original submission source)
"""

import os
import re
import sys
from collections import Counter
from difflib import SequenceMatcher

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
NEW = os.path.join(ROOT, 'paper', 'paper.tex')
DEFAULT_ORIG = r'..\programs\Resoumission\Resubmitted\Paper.tex'
# tolerated differences of the free edition (applied to the ORIGINAL side):
ORIG_STRIP = [
    r'\\footnote\{See supplementary material at \[URL will be inserted by AIP\] for \[SuppPub[123]\.zip[^}]*\}',
    r'\\docsection\{Figure Captions\}.*?\\end\{document\}',
]
# same, applied to the FREE side:
NEW_STRIP = [
    r'\\footnote\{The Matlab implementation.*?\}',
    r'\\footnote\{Provided in the companion repository.*?\}',
    r'\\docsection\{Multimedia and supplementary material\}.*?\\end\{itemize\}',
]
# purely cosmetic tokens introduced by the free template:
COSMETIC = {'\\raggedright', '\\small', '\\centering', '\\newpage',
            '\\setlength{\\tabcolsep}{3pt}', '\\end{document}',
            '\\bibliographystyle{plainnat}', '\\docsection{Tables}',
            '\\onecolumn', '\\hfill', '\\noindent', '\\columnbreak',
            '\\begin{minipage}[t]{0.485\\textwidth}', '\\begin{minipage}[t]{\\columnwidth}',
            '\\end{minipage}', '\\vspace{6pt}', '\\vspace{10pt}',
            '\\setlength{\\dbltextfloatsep}{10pt plus 2pt minus 2pt}',
            '\\renewcommand{\\arraystretch}{1.7}', '\\raggedbottom',
            '\\noindent\\begin{minipage}[t]{0.485\\textwidth}',
            '\\noindent\\begin{minipage}[t]{\\columnwidth}'}


def normalize(path, strip):
    s = open(path, encoding='utf-8', errors='replace').read().replace('\r\n', '\n')
    s = re.sub(r'(?<!\\)%.*', '', s)
    # footnotes : remplacees par des liens depot dans l'edition libre ->
    # on les retire des deux cotes pour ne comparer que le corps de texte
    s = re.sub(r'\\footnote\{(?:[^{}]|\{[^{}]*\})*\}', '', s)
    # legendes de tableaux en \captionof (grille 2x2) : ramenees a \caption
    s = s.replace('\\captionof{table}{', '\\caption{')
    # figure 3 legerement reduite (12cm -> 11.7cm) pour l'equation 8 en bas
    # de la colonne droite de sa page : neutralise pour la comparaison
    s = s.replace('width=11.7cm,height=11.7cm', 'width=12cm,height=12cm')
    # groupe d'espacement d'affichage autour de l'equation 8
    s = re.sub(r'\\begingroup\s*\\setlength\{\\abovedisplayskip\}\{[^}]*\}%?\s*'
               r'\\setlength\{\\belowdisplayskip\}\{[^}]*\}%?', '', s)
    s = s.replace('\\begingroup', '').replace('\\endgroup', '')
    # enveloppes de flottants (l'edition libre compose les tableaux en
    # flux deterministe) : on retire les wrappers des deux cotes
    s = re.sub(r'\\begin\{table\}(?:\[[!htbp]*\])?', '', s)
    s = s.replace('\\end{table}', '')
    s = s.replace('\\begin{center}', '').replace('\\end{center}', '')
    for pat in strip:
        s = re.sub(pat, '', s, flags=re.S)
    return re.sub(r'\s+', ' ', s)


def body(s):
    i = s.find(r'\section{\label{sec:1}')
    return s[i:] if i >= 0 else s


def check(label, items_o, items_n):
    missing = [x for x in items_o if x not in items_n]
    extra = [x for x in items_n if x not in items_o]
    ok = not missing and not extra
    print('[%s] %-34s origine=%d libre=%d' % ('OK  ' if ok else 'FAIL', label, len(items_o), len(items_n)))
    for x in missing[:5]:
        print('       absent de la version libre :', x[:130])
    for x in extra[:5]:
        print('       nouveau dans la version libre :', x[:130])
    return ok


def main(orig_path):
    orig = normalize(orig_path, ORIG_STRIP)
    new = normalize(NEW, NEW_STRIP)
    ok = True

    # correctif cosmetique documente : le tableau III d'origine declare
    # {ccccc} (5 colonnes) pour 4 colonnees utilisees ; l'edition libre
    # utilise {cccc}.  Seule la DERNIERE occurrence (tableau III) est
    # corrigee -- le tableau II utilise correctement ses 5 colonnes.
    quirk = r'\begin{tabular}{ccccc}'
    i = orig.rfind(quirk)
    orig = orig[:i] + quirk.replace('{ccccc}', '{cccc}') + orig[i + len(quirk):]

    # corps bruts (sans strips) pour les verifications ciblees
    raw_o = re.sub(r'(?<!\\)%.*', '', open(orig_path, encoding='utf-8',
                                           errors='replace').read())
    raw_n = re.sub(r'(?<!\\)%.*', '', open(NEW, encoding='utf-8',
                                           errors='replace').read())

    # ---- 1. word-level diff (cosmetic tokens tolerated) -------------------
    wo = [w for w in body(orig).split(' ') if w and w not in COSMETIC]
    wn = [w for w in body(new).split(' ') if w and w not in COSMETIC]
    sm = SequenceMatcher(None, wo, wn, autojunk=False)
    print('corps du texte : %d / %d mots, similarite %.4f' % (len(wo), len(wn), sm.ratio()))
    real_diffs = []
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == 'equal':
            continue
        a = ' '.join(wo[i1:i2])
        b = ' '.join(wn[j1:j2])
        if a in COSMETIC or b in COSMETIC or (not a and b in COSMETIC):
            continue
        # footnotes: origine = texte AIP supprime ; libre = lien depot ajoute
        if ('AIP' in a and 'repo' in b) or ('AIP' not in a and not a and 'repo' in b):
            continue
        if 'multimedia' in a.lower() or 'docsection' in a or 'bibliography' in a + b:
            continue
        real_diffs.append((a[:150], b[:150]))
    if real_diffs:
        ok = False
        print('[FAIL] %d difference(s) de contenu inattendue(s) :' % len(real_diffs))
        for a, b in real_diffs[:10]:
            print('   orig :', a)
            print('   libr :', b)
    else:
        print('[OK  ] aucune difference de contenu inattendue')

    # ---- 2. equations ------------------------------------------------------
    def equations(s):
        out = []
        for env in ('equation', 'align', 'subequations', 'empheq'):
            out += re.findall(r'\\begin\{%s\}.*?\\end\{%s\}' % (env, env), s)
        return out
    ok &= check('environnements d\'equations', equations(orig), equations(new))

    # ---- 3. table rows -----------------------------------------------------
    def rows(s):
        out = []
        for m in re.finditer(r'\\begin\{tabular\}\{[^}]*\}(.*?)\\end\{tabular\}', s, re.S):
            for line in m.group(1).split('\\\\'):
                r = re.sub(r'\s+', ' ', re.sub(r'\\hline', '', line)).strip(' &')
                if r:
                    out.append(r)
        return out
    ok &= check('lignes de tableaux', rows(orig), rows(new))

    # ---- 4. captions -------------------------------------------------------
    def captions(s):
        return re.findall(r'\\caption\{(?:\\label\{[^}]*\})?(.*?)\}\s*(?:\\raggedright|\\end\{figure)', s, re.S)
    ok &= check('legendes de figures', captions(orig), captions(new))

    # ---- 5. multimedia descriptions ----------------------------------------
    bo = raw_o[raw_o.find(r'\begin{document}'):]
    bn = raw_n[raw_n.find(r'\begin{document}'):]
    mo = re.findall(r'\\multimedia\{[^}]*\}\{(.*?)\}\\label', bo, re.S)
    mn = re.findall(r'\\multimedia\{[^}]*\}\{(.*?)\}\\label', bn, re.S)
    mo = [re.sub(r'\s+', ' ', m).strip() for m in mo]
    mn = [re.sub(r'\s+', ' ', m).strip() for m in mn]
    ok &= check('descriptions Mm.', mo, mn)

    # ---- 6. phonetic strings ------------------------------------------------
    ipa_o = Counter(re.findall(r'\\textipa\{[^}]*\}', bo))
    ipa_n = Counter(re.findall(r'\\textipa\{[^}]*\}', bn))
    # la version libre ne duplique pas les captions (liste de fin de soumission)
    for k, v in (ipa_o - ipa_n).items():
        ipa_o[k] = min(ipa_o[k], v + Counter(re.findall(re.escape(k),
                ' '.join(captions(orig))) and [1] or [])[k])
    diff = ipa_o - ipa_n
    if diff:
        ok = False
        print('[FAIL] \\textipa manquants :', dict(diff))
    else:
        print('[OK  ] caracteres phonetiques complets (%d occurrences)' % sum(ipa_n.values()))

    print('VERDICT :', 'TOUTES VERIFICATIONS PASSENT' if ok else 'ECHEC')
    return 0 if ok else 1


if __name__ == '__main__':
    orig = sys.argv[1] if len(sys.argv) > 1 else os.path.normpath(
        os.path.join(ROOT, DEFAULT_ORIG))
    if not os.path.exists(orig):
        sys.exit('Source original introuvable : %s' % orig)
    sys.exit(main(orig))
