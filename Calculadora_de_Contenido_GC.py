import streamlit as st
import re
from io import StringIO

def calcular_gc(secuencia):
    """
    Calcula el porcentaje de contenido GC en una secuencia de ADN.
    
    Args:
        secuencia (str): Secuencia de nucleótidos (ADN)
    
    Returns:
        float: Porcentaje de GC
    """
    secuencia = secuencia.upper().strip()
    # Eliminar caracteres no válidos (solo mantener ATCG)
    secuencia = re.sub(r'[^ATCG]', '', secuencia)
    
    total = len(secuencia)
    
    if total == 0:
        return 0.0, 0, 0, 0, 0
    
    count_g = secuencia.count('G')
    count_c = secuencia.count('C')
    count_a = secuencia.count('A')
    count_t = secuencia.count('T')
    
    count_gc = count_g + count_c
    porcentaje = (count_gc / total) * 100 if total > 0 else 0.0
    
    return round(porcentaje, 2), count_g, count_c, count_a, count_t

def parse_fasta(contenido):
    """Parsea el contenido de un archivo FASTA"""
    secuencias = {}
    current_id = ""
    current_seq = []
    
    for linea in contenido.splitlines():
        if linea.startswith('>'):
            if current_id:
                secuencias[current_id] = ''.join(current_seq)
            current_id = linea[1:].strip()
            current_seq = []
        else:
            current_seq.append(linea.strip())
    
    if current_id:
        secuencias[current_id] = ''.join(current_seq)
    
    return secuencias

# Configuración de la página
st.set_page_config(
    page_title="Calculadora de GC",
    page_icon="🧬",
    layout="centered",
)

# Título principal
st.title("🧬 Calculadora de Contenido GC")
st.caption("Herramienta para análisis de secuencias de ADN")

# Selección de modo de entrada
modo = st.radio("Seleccione el modo de entrada:", 
                ["Secuencia directa", "Archivo FASTA"], 
                horizontal=True)

resultados = []

if modo == "Secuencia directa":
    secuencia = st.text_area("Ingrese la secuencia de ADN:", 
                            height=200,
                            placeholder="Ej: ATGCGATACCTAG...")
    
    if st.button("Calcular GC"):
        if secuencia.strip():
            gc, g, c, a, t = calcular_gc(secuencia)
            resultados.append({
                "ID": "Secuencia ingresada",
                "GC": gc,
                "Longitud": len(re.sub(r'[^ATCG]', '', secuencia.upper())),
                "G": g,
                "C": c,
                "A": a,
                "T": t
            })
        else:
            st.warning("Por favor ingrese una secuencia de ADN")

else:
    archivo = st.file_uploader("Subir archivo FASTA", type=["fasta", "txt", "fa"])
    
    if archivo is not None:
        contenido = archivo.read().decode("utf-8")
        secuencias = parse_fasta(contenido)
        
        if not secuencias:
            st.error("No se encontraron secuencias válidas en el archivo")
        else:
            for seq_id, secuencia in secuencias.items():
                gc, g, c, a, t = calcular_gc(secuencia)
                resultados.append({
                    "ID": seq_id[:50] + "..." if len(seq_id) > 50 else seq_id,
                    "GC": gc,
                    "Longitud": len(secuencia),
                    "G": g,
                    "C": c,
                    "A": a,
                    "T": t
                })
            
            st.success(f"Se analizaron {len(secuencias)} secuencias")

# Mostrar resultados
if resultados:
    st.divider()
    st.subheader("Resultados del análisis")
    
    # Mostrar tabla resumen
    st.dataframe(
        resultados,
        column_config={
            "ID": "Identificador",
            "GC": st.column_config.NumberColumn(
                "% GC",
                format="%.2f %%",
                help="Porcentaje de Guanina + Citosina"
            ),
            "Longitud": "Longitud",
            "G": "Guanina",
            "C": "Citosina",
            "A": "Adenina",
            "T": "Timina"
        },
        hide_index=True,
        use_container_width=True
    )
    
    # Mostrar detalles para la primera secuencia
    if modo == "Secuencia directa":
        seq_data = resultados[0]
        col1, col2, col3 = st.columns(3)
        col1.metric("**% GC**", f"{seq_data['GC']}%")
        col2.metric("**Longitud**", seq_data['Longitud'])
        col3.metric("**Bases GC**", seq_data['G'] + seq_data['C'])
        
        # Gráfico de composición
        composicion = {
            'A': seq_data['A'],
            'T': seq_data['T'],
            'G': seq_data['G'],
            'C': seq_data['C']
        }
        
        st.bar_chart(composicion, color="#4CAF50")
        
        # Mostrar secuencia limpia
        st.caption("Secuencia procesada (solo bases ATCG):")
        st.code(secuencia[:500] + ("..." if len(secuencia) > 500 else ""))

# Información adicional
with st.expander("ℹ️ Acerca de esta herramienta"):
    st.markdown("""
    **¿Qué es el contenido de GC?**  
    El contenido de GC es el porcentaje de bases nitrogenadas en una molécula de ADN 
    que son Guanina (G) o Citosina (C). Es un indicador importante en biología molecular
    porque:
    - Secuencias con alto contenido GC son más estables térmicamente
    - Se correlaciona con regiones génicas importantes
    - Varía entre diferentes organismos y tipos de secuencias
    
    **Características de la calculadora:**
    - Elimina automáticamente caracteres no válidos
    - Soporta entrada directa o archivos FASTA
    - Proporciona análisis estadístico completo
    - Visualización gráfica de la composición de bases
    
    **Ejemplos de secuencias para probar:**
    ```
    Human BRCA1: GGATCCGATCGATCGATCGATAGCTAGCTAGCATCGATCGATCG
    Mitochondrial: ATGACTAGCTACGTACGTACGTACGTACGTACGATCGATCGTAC
    Bacterial: CGCGCGCGATATATATAGCGCGCGCTATATATAGCGCGCGC
    ```
    """)

# Pie de página
st.divider()
st.caption("Herramienta bioinformática desarrollada con Streamlit 🧬 | © 2023")
