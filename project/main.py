import pandas as pd

# 1. VERİYİ YÜKLE
# Not: Dosya yolunun doğruluğundan emin ol.
df = pd.read_csv(
    "data/pubchem_brd4_bioactivity_protein.csv",
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
    'Activity_Qualifier', 
    'Activity_Value', 
    'Compound_CID', 
    'Protein_Accession'
]

# Sadece ihtiyacımız olan sütunları al ve sayısal değeri/CID'si olmayanları sil
valid_metrics = ['IC50', 'Ki', 'Kd']
df_filtered = df[df['Activity_Type'].isin(valid_metrics)][selected_columns].copy()
df_filtered = df_filtered.dropna(subset=['Activity_Value', 'Compound_CID'])

# 3. ETİKETLEME FONKSİYONU (THRESHOLDING)
# Belirlediğin bilimsel kriterlere göre aktif/inaktif sınıflandırması yapar [cite: 17, 30]
# Değerlerin mikromolar olduğu varsayıldı.
def assign_label(row):
    m_type = row['Activity_Type']
    val = row['Activity_Value']
    
    # IC50 kriteri: < 1 uM Aktif, > 10 uM İnaktif
    if m_type == 'IC50':
        if val < 1.0: return 1
        if val > 10.0: return 0
    
    # Kd kriteri: < 500 nM (0.5 uM) Aktif, > 10 uM İnaktif
    elif m_type == 'Kd':
        if val < 0.5: return 1
        if val > 10.0: return 0
    
    # Ki kriteri: < 300 nM (0.3 uM) Aktif, > 10 uM İnaktif
    elif m_type == 'Ki':
        if val < 0.3: return 1
        if val > 10.0: return 0
    
    return None # Gri bölgede kalanlar

# Fonksiyonu uygula ve yeni Label sütununu oluştur. Yani etiketlemeyi uyguluyoruz
df_filtered['Label'] = df_filtered.apply(assign_label, axis=1)

# 5. Gri bölgeyi (None olanları) temizle ve duplikaları (mükerrer CID) sil
# Aynı bileşik hem aktif hem inaktif görünüyorsa en aktif olanı (1) tutarız
final_df = df_filtered.dropna(subset=['Label']).sort_values(by='Label', ascending=False)
final_df = final_df.drop_duplicates(subset=['Compound_CID'])

# Sonuçları ekrana yazdır
print("--- İşlem Tamamlandı ---")
print(f"Toplam Aktif (1) Sayısı: {len(final_df[final_df['Label'] == 1])}")
print(f"Toplam İnaktif (0) Sayısı: {len(final_df[final_df['Label'] == 0])}")

# Sonuçları kontrol et
print("--- TÜM VERİ SETİ SONUÇLARI ---")
print(final_df['Label'].value_counts())
print(f"Toplam benzersiz ve etiketli bileşik sayısı: {len(final_df)}")

# Temizlenmiş veriyi yeni bir CSV olarak kaydet
final_df.to_csv('brd4_egitim_verisi.csv', index=False)




