from sqlalchemy import create_engine, Column, Integer, String
from sqlalchemy.orm import declarative_base
from sqlalchemy.orm import relationship

db_path = "/n/core/Bioinformatics/analysis/CompBio/boweny/nf-Pipeline/Scundo3_v5/db/Secundo_Lite.db"
engine = create_engine(f"sqlite:///{db_path}", echo=True, future=True)

Base = declarative_base()

class order(Base):
    __tablename__ = "secundo_order"
    id = Column(Integer, primary_key=True)
    fcid = Column(String(50))
    molngid = Column(String(50))

    secundo_log = relationship("Secundo_Log", )