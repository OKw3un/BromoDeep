import pandas as pd
from pathlib import Path  # Modern ve temiz yol kontrolü için
import numpy as np

# ÇALIŞMA DİZİNİ KONTROLÜ (if-else bloğuna kadar olan 1-2-3 maddeleri için): 
# Path(__file__) kullanarak dosya yollarını çalıştırılan terminal dizinine göre değil, 
# kodun bilgisayardaki fiziksel konumuna göre belirliyoruz. 
# Böylece "cd project" yapsak bile data klasörü her zaman doğru bulunur.

# 1. KODUN OLDUĞU KLASÖRÜ BUL (Örn: .../bromodomain/project)
SCRIPT_DIR = Path(__file__).resolve().parent

# 2. BİR ÜST KLASÖRE ÇIK (Örn: .../bromodomain)
# Eğer data klasörü project'in içinde değil, yanındaysa .parent kullanmalıyız
BASE_DIR = SCRIPT_DIR.parent

# 3. YOLLARI BU ANA DİZİNE GÖRE AYARLA
# Bu sayede terminalde hangi klasörde olursan ol hata almazsın.
input_file = BASE_DIR / "data" / "pubchem_brd4_bioactivity_protein.csv"
output_file = BASE_DIR / "data" / "brd4_egitim_verisi.csv"

# Eğer dosya zaten varsa, işlemleri atla ve mevcut dosyayı yükle
if output_file.exists():
    print(f"--- '{output_file.name}' zaten mevcut. İşlemler atlanıyor. ---")
    final_df = pd.read_csv(output_file)
else:
    print(f"--- '{input_file.name}' bulunamadı. Veri işleme başlatılıyor... ---")

    # 1. VERİYİ YÜKLE
    # Not: Dosya yolunun doğruluğundan emin ol.
    df = pd.read_csv(
        input_file,
        sep=",",   
        engine="python",    #c parser yerine python parser kullanmayı seçtik. (python parser, c parser'a göre daha toleranslı ama yavaş.)
        on_bad_lines="skip", #bozuk bir satır görünce hata vermek yerine satır atlar
        encoding="utf-8" #Türkçe + özel karakterler için 
    )

    # 2. GEREKLİ SÜTUNLARI SEÇ VE TEMEL TEMİZLİK
    # Proje için kritik olan sütunları filtreliyoruz [cite: 14, 26]
    selected_columns = [
        'Activity',
        'Activity_Type',  
        'Activity_Value', 
        'Compound_CID', 
    ]

    # Sadece ihtiyacımız olan sütunları al ve sayısal değeri/CID'si olmayanları sil
    df_filtered = df[selected_columns].copy()
    # sadece IC50/Ki/Kd olanlar + Unspecified kalır
    df_filtered = df_filtered[
        df_filtered["Activity_Type"].isin(["IC50", "Ki", "Kd", "Unspecified"])
    ]
    # Activity_Value null ise sadece Unspecified ise sil
    df_filtered = df_filtered[
        ~(
            df_filtered["Activity_Value"].isna() &
            (df_filtered["Activity_Type"] == "Unspecified")
        )
    ]

    # 3. µM normalize (eğer veri zaten µM ise direkt kullan)
    df_filtered["uM"] = df_filtered["Activity_Value"]

    # 4. p-value dönüşümü
    df_filtered["p_value"] = 6 - np.log10(df_filtered["Activity_Value"])

    # 5. Label (1 µM threshold)
    df_filtered["Label"] = (df_filtered["uM"] <= 1.0).astype(int)
    df_filtered["Final_Activity"] = df_filtered["Activity"]
    mask = df_filtered["Activity"] == "Unspecified"
    df_filtered.loc[mask, "Final_Activity"] = np.where(
        df_filtered.loc[mask, "uM"] <= 1.0,
        "Active",
        "Inactive"
    )

    # duplicate cleanup + gereksiz sütunları kaldır
    final_df = (
    df_filtered
    .sort_values(by='p_value', ascending=False)
    .drop_duplicates(subset=['Compound_CID'])
    .drop(columns=["Label", "uM"])
    )

    # Kaydederken tam yolu kullan
    output_file.parent.mkdir(parents=True, exist_ok=True)
    final_df.to_csv(output_file, index=False)
    print("--- İşlem başarıyla tamamlandı ---")

print("İşlem tamamlandı.")

