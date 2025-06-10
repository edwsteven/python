# app.py
import streamlit as st
import gzip
import io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from collections import Counter, defaultdict
import base64
import os
import tempfile
import plotly.express as px

# Configuración de la página
st.set_page_config(
    page_title="FASTQ Quality Analyzer",
    page_icon=":dna:",
    layout="wide"
)

# Título de la aplicación
st.title('🧬 Analizador de Calidad FASTQ')
st.markdown("""
**Esta aplicación realiza análisis de calidad completo sobre archivos FASTQ/Fastq.gz incluyendo:**
- Estadísticas de calidad por posición
- Distribución de longitudes de lectura
- Conteo de bases ambiguas (N)
- Frecuencia de nucleótidos
- Detección de adaptadores
- Identificación de lecturas duplicadas
""")

# Sidebar para cargar archivo
with st.sidebar:
    st.header("Cargar archivo FASTQ")
    uploaded_file = st.file_uploader(
        "Suba su archivo FASTQ o FASTQ.GZ",
        type=['fastq', 'fq', 'fastq.gz', 'fq.gz']
    )
    
    st.markdown("---")
    st.info("""
    **Instrucciones:**
    1. Suba un archivo FASTQ (.fastq) o comprimido (.fastq.gz)
    2. Espere a que se complete el análisis (≈30 segundos)
    3. Explore los resultados y descargue el informe
    """)

# Función para analizar FASTQ
@st.cache_data(show_spinner="Analizando archivo FASTQ...")
def parse_fastq(uploaded_file):
    """Analiza el archivo FASTQ y extrae métricas básicas"""
    # Determinar si está comprimido
    is_gzipped = uploaded_file.name.endswith('.gz')
    
    # Crear un objeto archivo temporal
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fastq.gz" if is_gzipped else ".fastq") as tmp_file:
        tmp_file.write(uploaded_file.getvalue())
        tmp_path = tmp_file.name
    
    lengths = []
    qualities = []
    ambiguous_bases = []
    all_sequences = []
    nucleotide_counts = defaultdict(lambda: [0]*150)
    
    total_reads = 0
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    with open_func(tmp_path, mode) as handle:
        for i, record in enumerate(SeqIO.parse(handle, 'fastq')):
            # Solo procesar las primeras 10,000 lecturas para eficiencia
            if i >= 10000:
                break
                
            seq_len = len(record)
            lengths.append(seq_len)
            qual_scores = record.letter_annotations['phred_quality']
            qualities.append(qual_scores)
            ambiguous_bases.append(str(record.seq).count('N'))
            all_sequences.append(str(record.seq))
            
            # Frecuencia de nucleótidos por posición
            for pos, base in enumerate(str(record.seq)):
                if pos < 150:
                    nucleotide_counts[base][pos] += 1
            total_reads = i + 1
    
    # Eliminar archivo temporal
    os.unlink(tmp_path)
    
    return {
        'lengths': lengths,
        'qualities': qualities,
        'ambiguous_bases': ambiguous_bases,
        'sequences': all_sequences,
        'nucleotide_counts': nucleotide_counts,
        'total_reads': total_reads
    }

# Funciones de visualización
def plot_quality_scores(qualities):
    """Genera gráfico de calidad por posición"""
    max_len = max(len(q) for q in qualities)
    padded_qual = [q + [np.nan]*(max_len - len(q)) for q in qualities]
    df = pd.DataFrame(padded_qual)
    
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=df, color='skyblue', fliersize=1)
    plt.title('Calidad por Posición', fontsize=16)
    plt.xlabel('Posición en la lectura', fontsize=12)
    plt.ylabel('Calidad (Phred Score)', fontsize=12)
    plt.axhline(y=20, color='r', linestyle='--', alpha=0.7)
    plt.grid(True, alpha=0.3)
    st.pyplot(plt.gcf())
    plt.clf()

def plot_length_distribution(lengths):
    """Visualiza distribución de longitudes"""
    plt.figure(figsize=(10, 6))
    ax = sns.histplot(lengths, bins=50, kde=True, color='teal')
    plt.title('Distribución de Longitudes de Lectura', fontsize=16)
    plt.xlabel('Longitud (bp)', fontsize=12)
    plt.ylabel('Frecuencia', fontsize=12)
    plt.grid(True, alpha=0.3)
    st.pyplot(plt.gcf())
    plt.clf()

def plot_ambiguous_bases(ambiguous_bases):
    """Muestra distribución de bases N"""
    plt.figure(figsize=(10, 6))
    ax = sns.histplot(ambiguous_bases, bins=50, color='salmon')
    plt.title('Distribución de Bases Ambiguas (N) por Lectura', fontsize=16)
    plt.xlabel('Número de bases N', fontsize=12)
    plt.ylabel('Número de lecturas', fontsize=12)
    plt.grid(True, alpha=0.3)
    st.pyplot(plt.gcf())
    plt.clf()

def plot_nucleotide_frequency(nucleotide_counts):
    """Grafica frecuencia de nucleótidos por posición"""
    bases = ['A', 'T', 'C', 'G', 'N']
    positions = range(150)
    plt.figure(figsize=(10, 6))
    
    for base in bases:
        if base in nucleotide_counts:
            plt.plot(positions, nucleotide_counts[base][:150], label=base, linewidth=2.5)
    
    plt.title('Frecuencia de Nucleótidos por Posición', fontsize=16)
    plt.xlabel('Posición en la lectura', fontsize=12)
    plt.ylabel('Frecuencia (absoluta)', fontsize=12)
    plt.legend(title='Base')
    plt.grid(True, alpha=0.3)
    st.pyplot(plt.gcf())
    plt.clf()

def detect_adapters(sequences, adapters=['AGATCGGAAGAG', 'CTGTCTCTTATA']):
    """Detecta adaptadores comunes"""
    adapter_counts = {adapter: 0 for adapter in adapters}
    
    for seq in sequences:
        for adapter in adapters:
            if adapter in seq:
                adapter_counts[adapter] += 1
                
    # Convertir a porcentaje
    total = len(sequences)
    for adapter in adapter_counts:
        adapter_counts[adapter] = (adapter_counts[adapter] / total) * 100
        
    return adapter_counts

def find_duplicates(sequences):
    """Identifica lecturas duplicadas"""
    sequence_counter = Counter(sequences)
    duplicates = {seq: count for seq, count in sequence_counter.items() if count > 1}
    duplicate_percentage = (len(duplicates) / len(sequences)) * 100
    return duplicate_percentage

# Generar informe HTML
def generate_html_report(fastq_data, adapters, duplicates, filename):
    """Genera informe HTML con todos los resultados"""
    # Convertir imágenes a base64
    def fig_to_base64(fig):
        buf = io.BytesIO()
        fig.savefig(buf, format='png', bbox_inches='tight')
        buf.seek(0)
        return base64.b64encode(buf.read()).decode('utf-8')
    
    # Generar figuras
    figs = []
    for func in [plot_quality_scores, plot_length_distribution, 
                plot_ambiguous_bases, plot_nucleotide_frequency]:
        plt.figure()
        func(fastq_data)
        figs.append(plt.gcf())
        plt.close()
    
    fig1_base64 = fig_to_base64(figs[0])
    fig2_base64 = fig_to_base64(figs[1])
    fig3_base64 = fig_to_base64(figs[2])
    fig4_base64 = fig_to_base64(figs[3])
    
    # Estadísticas básicas
    stats_html = f"""
    <div class="stats-box">
        <h3>Estadísticas Generales</h3>
        <ul>
            <li>Total de lecturas: {fastq_data['total_reads']:,}</li>
            <li>Longitud promedio: {np.mean(fastq_data['lengths']):.1f} bp</li>
            <li>Calidad promedio: {np.mean([np.mean(q) for q in fastq_data['qualities']]):.1f}</li>
            <li>Bases ambiguas (N) promedio: {np.mean(fastq_data['ambiguous_bases']):.2f}</li>
            <li>Lecturas duplicadas: {duplicates:.2f}%</li>
        </ul>
    </div>
    """
    
    # Tabla de adaptadores
    adapters_html = "<table><tr><th>Adaptador</th><th>% de Detección</th></tr>"
    for adapter, percentage in adapters.items():
        adapters_html += f"<tr><td>{adapter}</td><td>{percentage:.2f}%</td></tr>"
    adapters_html += "</table>"
    
    # Construir HTML
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Reporte de Calidad FASTQ - {filename}</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; }}
            h1, h2, h3 {{ color: #2c3e50; }}
            .section {{ 
                background-color: #f8f9fa; 
                border-radius: 10px; 
                padding: 20px; 
                margin-bottom: 30px;
                box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            }}
            .grid-container {{
                display: grid;
                grid-template-columns: repeat(2, 1fr);
                gap: 20px;
            }}
            .grid-item {{ background: white; padding: 15px; border-radius: 5px; }}
            .stats-box {{ 
                background: #e8f4f8; 
                padding: 15px; 
                border-left: 4px solid #3498db;
                margin-bottom: 20px;
            }}
            table {{ 
                width: 100%; 
                border-collapse: collapse; 
                margin: 15px 0;
            }}
            th, td {{ 
                padding: 10px; 
                text-align: left; 
                border-bottom: 1px solid #ddd;
            }}
            img {{ max-width: 100%; height: auto; border: 1px solid #ddd; }}
        </style>
    </head>
    <body>
        <h1>Reporte de Calidad FASTQ</h1>
        <h2>Archivo: {filename}</h2>
        
        <div class="section">
            <h3>Resumen de Calidad</h3>
            {stats_html}
        </div>
        
        <div class="section">
            <h3>Detección de Adaptadores</h3>
            {adapters_html}
        </div>
        
        <div class="section">
            <h3>Visualizaciones</h3>
            <div class="grid-container">
                <div class="grid-item">
                    <h4>Calidad por Posición</h4>
                    <img src="data:image/png;base64,{fig1_base64}" alt="Quality Scores">
                </div>
                <div class="grid-item">
                    <h4>Distribución de Longitudes</h4>
                    <img src="data:image/png;base64,{fig2_base64}" alt="Length Distribution">
                </div>
                <div class="grid-item">
                    <h4>Bases Ambiguas (N)</h4>
                    <img src="data:image/png;base64,{fig3_base64}" alt="Ambiguous Bases">
                </div>
                <div class="grid-item">
                    <h4>Frecuencia de Nucleótidos</h4>
                    <img src="data:image/png;base64,{fig4_base64}" alt="Nucleotide Frequency">
                </div>
            </div>
        </div>
        
        <div class="section">
            <h3>Métricas Detalladas</h3>
            <p><strong>Calidad Promedio:</strong> {np.mean([np.mean(q) for q in fastq_data['qualities']]):.1f}</p>
            <p><strong>Longitud Mínima:</strong> {np.min(fastq_data['lengths'])} bp</p>
            <p><strong>Longitud Máxima:</strong> {np.max(fastq_data['lengths'])} bp</p>
            <p><strong>% Lecturas con N:</strong> {np.mean([1 if n > 0 else 0 for n in fastq_data['ambiguous_bases']]) * 100:.2f}%</p>
        </div>
    </body>
    </html>
    """
    
    return html_content

# Mostrar contenido principal
if uploaded_file is not None:
    # Analizar archivo FASTQ
    with st.spinner('Analizando archivo FASTQ...'):
        fastq_data = parse_fastq(uploaded_file)
    
    # Mostrar estadísticas básicas
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Lecturas totales", f"{fastq_data['total_reads']:,}")
    col2.metric("Longitud promedio", f"{np.mean(fastq_data['lengths']):.1f} bp")
    col3.metric("Calidad promedio", f"{np.mean([np.mean(q) for q in fastq_data['qualities']]):.1f}")
    col4.metric("Bases N promedio", f"{np.mean(fastq_data['ambiguous_bases']):.2f}")
    
    st.divider()
    
    # Análisis adicional
    adapters = detect_adapters(fastq_data['sequences'])
    duplicates = find_duplicates(fastq_data['sequences'])
    
    # Mostrar resultados en pestañas
    tab1, tab2, tab3, tab4 = st.tabs([
        "Calidad por Posición", 
        "Longitudes", 
        "Bases Ambiguas", 
        "Frecuencia Nucleótidos"
    ])
    
    with tab1:
        st.subheader("Distribución de Calidad por Posición")
        plot_quality_scores(fastq_data['qualities'])
        
    with tab2:
        st.subheader("Distribución de Longitudes de Lectura")
        plot_length_distribution(fastq_data['lengths'])
        
    with tab3:
        st.subheader("Distribución de Bases Ambiguas (N)")
        plot_ambiguous_bases(fastq_data['ambiguous_bases'])
        
    with tab4:
        st.subheader("Frecuencia de Nucleótidos por Posición")
        plot_nucleotide_frequency(fastq_data['nucleotide_counts'])
    
    st.divider()
    
    # Resultados adicionales
    st.subheader("Resultados Adicionales")
    
    col5, col6 = st.columns(2)
    
    with col5:
        st.markdown("**Detección de Adaptadores**")
        for adapter, percentage in adapters.items():
            st.progress(float(percentage/100), text=f"{adapter}: {percentage:.2f}%")
    
    with col6:
        st.markdown("**Lecturas Duplicadas**")
        st.metric("Porcentaje de duplicados", f"{duplicates:.2f}%")
    
    st.divider()
    
    # Generar y descargar informe
    st.subheader("Reporte Completo")
    
    if st.button("Generar Informe HTML"):
        with st.spinner("Generando reporte..."):
            report = generate_html_report(
                fastq_data, 
                adapters, 
                duplicates, 
                uploaded_file.name
            )
            
            st.download_button(
                label="Descargar Reporte",
                data=report,
                file_name="fastq_report.html",
                mime="text/html"
            )
else:
    st.info("Por favor, suba un archivo FASTQ para comenzar el análisis")
    st.image("https://upload.wikimedia.org/wikipedia/commons/thumb/0/0d/FASTQ_format_Illumina.png/800px-FASTQ_format_Illumina.png", 
             caption="Formato FASTQ típico", use_column_width=True)

# Notas al pie
st.divider()
st.caption("FASTQ Quality Analyzer v1.0 | Desarrollado con Streamlit y Biopython")