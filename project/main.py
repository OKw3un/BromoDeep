import pandas as pd

df = pd.read_csv(
    "data/Kitap1.csv",
    sep=",",   
    engine="python",    #c parser yerine python parser kullanmayı seçtik. (python parser, c parser'a göre daha toleranslı ama yavaş.)
    on_bad_lines="skip", #bozuk bir satır görünce hata vermek yerine satır atlar
    encoding="utf-8" #Türkçe + özel karakterler için 
)

"""
pd.set_option("display.max_columns", None)
pd.set_option("display.width", None)
print(df)
print(df.shape)
"""

# 2. Sadece senin istediğin 5 sütunu seçelim
selected_columns = [
    'Activity',
    'Activity_Type', 
    'Activity_Qualifier', 
    'Activity_Value', 
    'Compound_CID', 
    'Protein_Accession'
]

# Yeni bir DataFrame oluştur
egitim_df = df[selected_columns]

# --- ANALİZ VE GÖSTERİM ---

# 3. Önce ilk 100 satırı filtrele, sonra içinden rastgele 10 tanesini seç
# .iloc[:100] -> ilk 100 satırı alır
# .sample(10) -> içinden rastgele 10 tanesini çeker
random_sample = df[selected_columns].iloc[:100].sample(10)
print("--- İlk 100 Satır İçinden Rastgele Seçilen 10 Örnek ---")
print(random_sample)

# Veri seti hakkında kısa özet (Kaç satır dolu, veri tipleri ne?)
print("\n--- Veri Seti Özeti ---")
print(egitim_df.info())

# Kaç tane boş (NaN) değer olduğunu göster
print("\n--- Eksik Değer Sayısı ---")
print(egitim_df.isnull().sum())

# Activity_Type dağılımını göster (Kaç IC50, kaç Tm var?)
print("\n--- Aktivite Türü Dağılımı ---")
print(egitim_df['Activity_Type'].value_counts())

