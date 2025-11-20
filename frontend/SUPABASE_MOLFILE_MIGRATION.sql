-- Migration: Add molfile field to molecules table
-- Run this in Supabase SQL Editor

-- Add molfile column (nullable, as existing molecules won't have it)
ALTER TABLE molecules 
ADD COLUMN IF NOT EXISTS molfile TEXT;

-- Add index for faster searches (optional, but recommended if you search by molfile)
-- CREATE INDEX IF NOT EXISTS idx_molecules_molfile ON molecules USING gin(to_tsvector('english', molfile));

-- Update existing molecules to generate molfiles (optional)
-- This would require calling the /convert/smiles endpoint for each molecule
-- You can do this programmatically or leave it for users to regenerate

-- Verify the column was added
SELECT column_name, data_type, is_nullable 
FROM information_schema.columns 
WHERE table_name = 'molecules' AND column_name = 'molfile';


