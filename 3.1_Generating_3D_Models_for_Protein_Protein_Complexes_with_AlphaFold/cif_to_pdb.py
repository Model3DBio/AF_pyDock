#!/usr/bin/env python3
import argparse
from Bio.PDB import MMCIFParser, PDBIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def _as_list(x):
    if x is None:
        return []
    return x if isinstance(x, list) else [x]


def _val_at(lst, i):
    if i >= len(lst):
        return ""
    v = lst[i]
    if v in (None, ".", "?", ""):
        return ""
    return str(v)


def _pick_per_row(label_list, auth_list, i):
    v = _val_at(label_list, i)
    return v if v else _val_at(auth_list, i)


def _norm_atom_name(x: str) -> str:
    return (x or "").strip()


def _norm_resname(x: str) -> str:
    return (x or "").strip()


def _parse_pdb_atom_index(pdb_lines):
    """
    (chain, resseq, icode, resname, atomname, altloc) -> serial (int)
    """
    idx = {}
    for line in pdb_lines:
        if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
            continue
        try:
            serial = int(line[6:11])
        except ValueError:
            continue

        atomname = line[12:16].strip()
        altloc = line[16].strip()
        resname = line[17:20].strip()
        chain = line[21].strip()
        resseq = line[22:26].strip()
        icode = line[26].strip()

        key = (chain, resseq, icode, resname, atomname, altloc)
        idx[key] = serial

        key_no_alt = (chain, resseq, icode, resname, atomname, "")
        idx.setdefault(key_no_alt, serial)

    return idx


def _get_struct_conn_rows(cif_dict):
    """
    Lee filas _struct_conn normalizadas (label por fila si no '.', si no auth).
    Claves ins_code/altloc según tu CIF: pdbx_ptnr*_...
    """
    conn_type = _as_list(cif_dict.get("_struct_conn.conn_type_id"))
    if not conn_type:
        return []

    # ptnr1
    p1_lab_asym = _as_list(cif_dict.get("_struct_conn.ptnr1_label_asym_id"))
    p1_auth_asym = _as_list(cif_dict.get("_struct_conn.ptnr1_auth_asym_id"))
    p1_lab_seq = _as_list(cif_dict.get("_struct_conn.ptnr1_label_seq_id"))
    p1_auth_seq = _as_list(cif_dict.get("_struct_conn.ptnr1_auth_seq_id"))
    p1_lab_comp = _as_list(cif_dict.get("_struct_conn.ptnr1_label_comp_id"))
    p1_auth_comp = _as_list(cif_dict.get("_struct_conn.ptnr1_auth_comp_id"))
    p1_lab_atom = _as_list(cif_dict.get("_struct_conn.ptnr1_label_atom_id"))
    p1_auth_atom = _as_list(cif_dict.get("_struct_conn.ptnr1_auth_atom_id"))

    # ptnr2
    p2_lab_asym = _as_list(cif_dict.get("_struct_conn.ptnr2_label_asym_id"))
    p2_auth_asym = _as_list(cif_dict.get("_struct_conn.ptnr2_auth_asym_id"))
    p2_lab_seq = _as_list(cif_dict.get("_struct_conn.ptnr2_label_seq_id"))
    p2_auth_seq = _as_list(cif_dict.get("_struct_conn.ptnr2_auth_seq_id"))
    p2_lab_comp = _as_list(cif_dict.get("_struct_conn.ptnr2_label_comp_id"))
    p2_auth_comp = _as_list(cif_dict.get("_struct_conn.ptnr2_auth_comp_id"))
    p2_lab_atom = _as_list(cif_dict.get("_struct_conn.ptnr2_label_atom_id"))
    p2_auth_atom = _as_list(cif_dict.get("_struct_conn.ptnr2_auth_atom_id"))

    # claves correctas (como en tu fragmento)
    p1_ins = _as_list(cif_dict.get("_struct_conn.pdbx_ptnr1_PDB_ins_code"))
    p2_ins = _as_list(cif_dict.get("_struct_conn.pdbx_ptnr2_PDB_ins_code"))
    p1_alt = _as_list(cif_dict.get("_struct_conn.pdbx_ptnr1_label_alt_id"))
    p2_alt = _as_list(cif_dict.get("_struct_conn.pdbx_ptnr2_label_alt_id"))

    rows = []
    for i in range(len(conn_type)):
        ctype = (conn_type[i] or "").lower()
        p1 = {
            "chain": _pick_per_row(p1_lab_asym, p1_auth_asym, i),
            "resseq": _pick_per_row(p1_lab_seq, p1_auth_seq, i),
            "icode": _val_at(p1_ins, i),
            "altloc": _val_at(p1_alt, i),
            "resname": _pick_per_row(p1_lab_comp, p1_auth_comp, i),
            "atomname": _pick_per_row(p1_lab_atom, p1_auth_atom, i),
        }
        p2 = {
            "chain": _pick_per_row(p2_lab_asym, p2_auth_asym, i),
            "resseq": _pick_per_row(p2_lab_seq, p2_auth_seq, i),
            "icode": _val_at(p2_ins, i),
            "altloc": _val_at(p2_alt, i),
            "resname": _pick_per_row(p2_lab_comp, p2_auth_comp, i),
            "atomname": _pick_per_row(p2_lab_atom, p2_auth_atom, i),
        }
        rows.append({"type": ctype, "p1": p1, "p2": p2})
    return rows


def _lookup_serial(atom_index, p):
    chain = (p.get("chain") or "").strip()
    resseq = (p.get("resseq") or "").strip()
    icode = (p.get("icode") or "").strip()
    resname = (p.get("resname") or "").strip()
    atomname = (p.get("atomname") or "").strip()
    altloc = (p.get("altloc") or "").strip()

    key = (chain, resseq, icode, resname, atomname, altloc)
    s = atom_index.get(key)

    if s is None and altloc:
        key2 = (chain, resseq, icode, resname, atomname, "")
        s = atom_index.get(key2)
    if s is not None:
        return s

    for (ch, rs, ic, rn, an, al), serial in atom_index.items():
        if ch == chain and rs == resseq and ic == icode and an == atomname:
            if altloc and al not in ("", altloc):
                continue
            return serial
    return None


def _format_link_record(p1, p2):
    atom1 = _norm_atom_name(p1["atomname"])
    res1  = _norm_resname(p1["resname"])
    ch1   = ((p1.get("chain") or " ").strip()[:1] or " ")
    rs1   = (p1.get("resseq") or "").strip()
    ic1   = ((p1.get("icode") or " ").strip()[:1] or " ")

    atom2 = _norm_atom_name(p2["atomname"])
    res2  = _norm_resname(p2["resname"])
    ch2   = ((p2.get("chain") or " ").strip()[:1] or " ")
    rs2   = (p2.get("resseq") or "").strip()
    ic2   = ((p2.get("icode") or " ").strip()[:1] or " ")

    line = [" "] * 80
    line[0:6] = list("LINK  ")

    def put(pos1, text):
        start = pos1 - 1
        for j, c in enumerate(text):
            if 0 <= start + j < 80:
                line[start + j] = c

    put(13, atom1.rjust(4))
    put(18, res1.rjust(3))
    put(22, ch1)
    put(23, rs1.rjust(4))
    put(27, ic1)

    put(43, atom2.rjust(4))
    put(48, res2.rjust(3))
    put(52, ch2)
    put(53, rs2.rjust(4))
    put(57, ic2)

    return "".join(line)


def _format_conect_record(serial_a, serial_b):
    return f"CONECT{serial_a:5d}{serial_b:5d}".ljust(80)


def _format_modres_record(idcode4, resname, chain, resseq, icode, stdres, comment):
    idc = (idcode4 or "1CIF")[:4].ljust(4)
    rn = (resname or "").strip()[:3].rjust(3)
    ch = ((chain or " ").strip()[:1] or " ")
    rs = (str(resseq) if resseq is not None else "").strip()
    ic = ((icode or " ").strip()[:1] or " ")
    sr = (stdres or "").strip()[:3].rjust(3)
    cm = (comment or "").strip()[:41]

    line = [" "] * 80
    line[0:6] = list("MODRES")

    def put(pos1, text):
        start = pos1 - 1
        for j, c in enumerate(text):
            if 0 <= start + j < 80:
                line[start + j] = c

    put(8, idc)
    put(13, rn)
    put(17, ch)
    put(19, rs.rjust(4))
    put(23, ic)
    put(25, sr)
    put(30, cm)

    return "".join(line)


def _format_ssbond_record(n, chain1, resseq1, icode1, chain2, resseq2, icode2):
    ch1 = ((chain1 or " ").strip()[:1] or " ")
    ch2 = ((chain2 or " ").strip()[:1] or " ")
    rs1 = (str(resseq1) if resseq1 is not None else "").strip()
    rs2 = (str(resseq2) if resseq2 is not None else "").strip()
    ic1 = ((icode1 or " ").strip()[:1] or " ")
    ic2 = ((icode2 or " ").strip()[:1] or " ")

    line = [" "] * 80
    line[0:6] = list("SSBOND")

    def put(pos1, text):
        start = pos1 - 1
        for j, c in enumerate(text):
            if 0 <= start + j < 80:
                line[start + j] = c

    put(8, f"{n:>3d}")
    put(12, "CYS")
    put(16, ch1)
    put(18, rs1.rjust(4))
    put(22, ic1)

    put(26, "CYS")
    put(30, ch2)
    put(32, rs2.rjust(4))
    put(36, ic2)

    return "".join(line)


def _dedupe(lines):
    seen = set()
    out = []
    for ln in lines:
        if ln not in seen:
            out.append(ln)
            seen.add(ln)
    return out


def _infer_ssbonds_from_distance(structure, max_dist=2.2):
    """
    Detecta puentes disulfuro por distancia SG-SG (Å).
    Devuelve pares (chain, resseq, icode) para escribir SSBOND.
    """
    sg_atoms = []  # (chain, resseq(str), icode(str), coord(np-like))
    for model in structure:
        for chain in model:
            for res in chain:
                if res.get_resname().strip() != "CYS":
                    continue
                if "SG" not in res:
                    continue
                atom = res["SG"]
                hetflag, resseq, icode = res.id[0], res.id[1], (res.id[2] or " ")
                # Normaliza icode a '' o letra
                ic = icode.strip()
                sg_atoms.append((chain.id, str(resseq), ic, atom.coord))

    pairs = []
    n = len(sg_atoms)
    for i in range(n):
        ch1, rs1, ic1, c1 = sg_atoms[i]
        for j in range(i + 1, n):
            ch2, rs2, ic2, c2 = sg_atoms[j]
            dx = float(c1[0] - c2[0])
            dy = float(c1[1] - c2[1])
            dz = float(c1[2] - c2[2])
            d = (dx * dx + dy * dy + dz * dz) ** 0.5
            if d <= max_dist:
                pairs.append(((ch1, rs1, ic1), (ch2, rs2, ic2)))
    return pairs


def convert_cif_to_pdb_with_records(
    cif_file,
    pdb_file,
    add_link=True,
    add_conect=True,
    add_modres=True,
    add_ssbond=True,
    ssbond_max_dist=2.2,
    idcode=None,
    auth_chains=True,
    auth_residues=True,
):
    # 1) mmCIF -> PDB
    parser = MMCIFParser(QUIET=True, auth_chains=auth_chains, auth_residues=auth_residues)
    structure = parser.get_structure("ID", cif_file)

    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)

    # 2) Leer PDB para índice
    with open(pdb_file, "r", encoding="utf-8") as f:
        pdb_lines = f.read().splitlines()
    atom_index = _parse_pdb_atom_index(pdb_lines)

    # 3) Leer CIF dict
    cif_dict = MMCIF2Dict(cif_file)
    if not idcode:
        ent = cif_dict.get("_entry.id")
        if isinstance(ent, list):
            ent = ent[0] if ent else None
        idcode = (str(ent).strip()[:4] if ent else "1CIF")

    rows = _get_struct_conn_rows(cif_dict)

    link_lines = []
    conect_lines = []
    modres_lines = []
    ssbond_pairs = []  # ((ch, rs, ic), (ch, rs, ic))
    glyco_sites = set()  # (ch, rs, ic)

    # 4) LINK/CONECT + MODRES desde struct_conn (si existe)
    for r in rows:
        ctype = (r["type"] or "").lower()
        p1, p2 = r["p1"], r["p2"]

        # LINK/CONECT: covale/disulf (si existiese disulf)
        if any(k in ctype for k in ("covale", "disulf")):
            s1 = _lookup_serial(atom_index, p1)
            s2 = _lookup_serial(atom_index, p2)
            if s1 is not None and s2 is not None:
                if add_link:
                    link_lines.append(_format_link_record(p1, p2))
                if add_conect:
                    conect_lines.append(_format_conect_record(s1, s2))
                    conect_lines.append(_format_conect_record(s2, s1))

        # MODRES: ASN ND2 -- NAG C1 (covalente)
        if add_modres and "covale" in ctype:
            rn1 = _norm_resname(p1.get("resname"))
            rn2 = _norm_resname(p2.get("resname"))
            a1 = _norm_atom_name(p1.get("atomname"))
            a2 = _norm_atom_name(p2.get("atomname"))

            is_asn_to_nag = (rn1 == "ASN" and a1 == "ND2" and rn2 == "NAG" and a2 == "C1")
            is_nag_to_asn = (rn2 == "ASN" and a2 == "ND2" and rn1 == "NAG" and a1 == "C1")

            if is_asn_to_nag:
                glyco_sites.add((p1.get("chain", ""), str(p1.get("resseq", "")), p1.get("icode", "") or ""))
            elif is_nag_to_asn:
                glyco_sites.add((p2.get("chain", ""), str(p2.get("resseq", "")), p2.get("icode", "") or ""))

    # 5) MODRES líneas
    if add_modres and glyco_sites:
        def _sort_key(t):
            ch, rs, ic = t
            try:
                rsi = int(rs)
            except Exception:
                rsi = 10**9
            return (ch, rsi, ic)

        for (ch, rs, ic) in sorted(glyco_sites, key=_sort_key):
            modres_lines.append(
                _format_modres_record(
                    idcode4=idcode,
                    resname="ASN",
                    chain=ch,
                    resseq=rs,
                    icode=ic,
                    stdres="ASN",
                    comment="GLYCOSYLATION SITE",
                )
            )

    # 6) SSBOND por distancia (independiente de struct_conn)
    if add_ssbond:
        inferred = _infer_ssbonds_from_distance(structure, max_dist=ssbond_max_dist)
        ssbond_pairs.extend(inferred)

    # dedupe pares SSBOND (sin importar orden)
    uniq = set()
    ordered = []
    for (a, b) in ssbond_pairs:
        key = tuple(sorted([a, b]))
        if key not in uniq:
            uniq.add(key)
            ordered.append((a, b))

    ssbond_lines = []
    if add_ssbond and ordered:
        for i, (a, b) in enumerate(ordered, start=1):
            (ch1, rs1, ic1) = a
            (ch2, rs2, ic2) = b
            ssbond_lines.append(_format_ssbond_record(i, ch1, rs1, ic1, ch2, rs2, ic2))

    # Dedupe final
    ssbond_lines = _dedupe(ssbond_lines)
    modres_lines = _dedupe(modres_lines)
    link_lines = _dedupe(link_lines)
    conect_lines = _dedupe(conect_lines)

    # 7) Insertar cabecera antes de ATOM/HETATM/MODEL
    insert_idx = 0
    for i, line in enumerate(pdb_lines):
        if line.startswith(("ATOM  ", "HETATM", "MODEL ")):
            insert_idx = i
            break
        insert_idx = i + 1

    # 8) Insertar CONECT antes de END
    end_idx = None
    for i in range(len(pdb_lines) - 1, -1, -1):
        if pdb_lines[i].startswith("END"):
            end_idx = i
            break
    if end_idx is None:
        end_idx = len(pdb_lines)

    new_lines = []
    new_lines.extend(pdb_lines[:insert_idx])

    # Orden típico
    if ssbond_lines:
        new_lines.extend(ssbond_lines)
    if modres_lines:
        new_lines.extend(modres_lines)
    if link_lines:
        new_lines.extend(link_lines)

    new_lines.extend(pdb_lines[insert_idx:end_idx])

    if conect_lines:
        new_lines.extend(conect_lines)

    new_lines.extend(pdb_lines[end_idx:])

    with open(pdb_file, "w", encoding="utf-8") as f:
        f.write("\n".join(new_lines) + "\n")


def main():
    ap = argparse.ArgumentParser(
        description="Convert CIF->PDB with LINK/CONECT from _struct_conn + MODRES (ASN-NAG) + SSBOND by SG-SG distance."
    )
    ap.add_argument("cif_file", help="Input mmCIF file (.cif)")
    ap.add_argument("pdb_file", help="Output PDB file (.pdb)")
    ap.add_argument("--no-link", action="store_true", help="Do not write LINK records")
    ap.add_argument("--no-conect", action="store_true", help="Do not write CONECT records")
    ap.add_argument("--no-modres", action="store_true", help="Do not write MODRES records")
    ap.add_argument("--no-ssbond", action="store_true", help="Do not write SSBOND records")
    ap.add_argument("--ssbond-max-dist", type=float, default=2.2, help="Max SG-SG distance (Å) to call disulfide (default: 2.2)")
    ap.add_argument("--idcode", type=str, default=None, help="Override PDB idCode (4 chars) used in MODRES")
    ap.add_argument("--no-auth", action="store_true", help="Do not force auth_chains/auth_residues in MMCIFParser")
    args = ap.parse_args()

    use_auth = not args.no_auth
    convert_cif_to_pdb_with_records(
        args.cif_file,
        args.pdb_file,
        add_link=not args.no_link,
        add_conect=not args.no_conect,
        add_modres=not args.no_modres,
        add_ssbond=not args.no_ssbond,
        ssbond_max_dist=args.ssbond_max_dist,
        idcode=args.idcode,
        auth_chains=use_auth,
        auth_residues=use_auth,
    )


if __name__ == "__main__":
    main()
